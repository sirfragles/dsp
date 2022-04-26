#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "matrix4.h"
#include "biquad.h"
#include "util.h"

#define ADAPT_TIME          300.0
#define RISE_TIME_FAST       33.0
#define RISE_TIME_SLOW      100.0
#define NORM_TIME           160.0
#define EVENT_THRESH          1.5
#define EVENT_END_THRESH      0.2
#define NORM_MAX              3.2
#define NORM_CROSSFEED        0.5
#define EVENT_SAMPLE_TIME    16.5
#define EVENT_MAX_HOLD_TIME 200.0
#define EVENT_MIN_HOLD_TIME  50.0

struct ewma_state {
	double c0, c1, m0;
};

enum event_state {
	EVENT_STATE_NONE = 0,
	EVENT_STATE_SAMPLE,
	EVENT_STATE_HOLD,
};

enum {
	EVENT_FLAG_L = 1<<0,
	EVENT_FLAG_R = 1<<1,
	EVENT_FLAG_USE_ORD = 1<<2,
	EVENT_FLAG_END = 1<<4,
};

struct ap1_state {
	sample_t c0;
	sample_t i0, o0;
};

struct ap2_state {
	sample_t c0, c1;
	sample_t i0, i1, o0, o1;
};

struct ap3_state {
	struct ap2_state ap2;
	struct ap1_state ap1;
};

struct cap5_state {
	struct ap2_state a1;
	struct ap3_state a2;
};

struct filter_bank {
	struct cap5_state b[5];
	struct ap2_state ap[6];
	sample_t s[6];
};

struct matrix4_band {
	enum event_state ev;
	int ev_flags;
	struct ewma_state r_fast[4], r_slow[2], adapt[4], norm[4], smooth[2], avg[2], drift[4];
	ssize_t t;
	ssize_t ord_count, diff_count, early_count;
	double lr, cs, lr_adapt, cs_adapt, fl_boost, fr_boost;
};

struct matrix4_state {
	int c0, c1, has_output, is_draining, disable, show_status, do_dir_boost;
	struct filter_bank fb[2];
	struct matrix4_band band[5];
	sample_t **bufs;
	sample_t *fb_buf[2][6];
	sample_t norm_mult, surr_mult;
	struct ewma_state band_d_sum_pwr_env[6];
	struct biquad_state band0_hp;  /* for computing band0 weight for directional boost */
	ssize_t len, p, drain_frames;
	ssize_t ev_sample_frames, ev_max_hold_frames, ev_min_hold_frames;
};

static void ewma_init(struct ewma_state *state, double fs, double f0)
{
	const double w0 = 2.0*M_PI*f0/fs;
	const double y = 1.0-cos(w0);
	const double a = -y+sqrt(y*y+2.0*y);
	state->c0 = a;
	state->c1 = 1.0-a;
	state->m0 = 0.0;
}

static double ewma_run(struct ewma_state *state, double s)
{
	const double r = state->c0*s + state->c1*state->m0;
	state->m0 = r;
	return r;
}

static double ewma_run_set_max(struct ewma_state *state, double s)
{
	if (s >= state->m0) s = ewma_run(state, s);
	else state->m0 = s;
	return s;
}

static double ewma_set(struct ewma_state *state, double s)
{
	state->m0 = s;
	return s;
}

static double ewma_get_last(struct ewma_state *state)
{
	return state->m0;
}

static sample_t ap1_run(struct ap1_state *state, sample_t s)
{
	sample_t r = state->c0 * (s - state->o0)
		+ state->i0;

	state->i0 = s;
	state->o0 = r;

	return r;
}

static sample_t ap2_run(struct ap2_state *state, sample_t s)
{
	sample_t r = state->c1 * (s - state->o1)
		+ state->c0 * (state->i0 - state->o0)
		+ state->i1;

	state->i1 = state->i0;
	state->i0 = s;

	state->o1 = state->o0;
	state->o0 = r;

	return r;
}

static sample_t ap3_run(struct ap3_state *state, sample_t s)
{
	return ap1_run(&state->ap1, ap2_run(&state->ap2, s));
}

/* doubly complementary 5th-order butterworth filters implemented as the
   sum (lowpass) and difference (highpass) of two allpass sections */
static void cap5_init(struct cap5_state *state, double fs, double fc)
{
	int i;
	double theta;
	double fc_w = fs/M_PI * tan(M_PI*fc/fs);  /* pre-warped corner frequency */
	double complex p[3];  /* first two have a complex conjugate (not stored), third is real */

	for (i = 0; i < 3; ++i) {
		theta = (2.0*(i+1)-1.0)*M_PI/10.0;
		p[i] = -sin(theta) + cos(theta)*I;  /* normalized pole in s-plane */
		p[i] = p[i]*2.0*M_PI*fc_w;          /* scale */
		p[i] = (1.0 + p[i]/(2.0*fs)) / (1.0 - p[i]/(2.0*fs));  /* bilinear transform */
		//LOG_FMT(LL_NORMAL, "fc=%gHz: p[%d] = %f%+fi", fc, i, creal(p[i]), cimag(p[i]));
	}

	state->a2.ap2.c0 = -2.0*creal(p[0]);
	state->a2.ap2.c1 = pow(creal(p[0]), 2.0) + pow(cimag(p[0]), 2.0);

	state->a1.c0 = -2.0*creal(p[1]);
	state->a1.c1 = pow(creal(p[1]), 2.0) + pow(cimag(p[1]), 2.0);

	state->a2.ap1.c0 = -creal(p[2]);

	//LOG_FMT(LL_NORMAL, "fc=%gHz: a1: c0=%g  c1=%g", fc, state->a1.c0, state->a1.c1);
	//LOG_FMT(LL_NORMAL, "fc=%gHz: a2.ap2: c0=%g  c1=%g", fc, state->a2.ap2.c0, state->a2.ap2.c1);
	//LOG_FMT(LL_NORMAL, "fc=%gHz: a2.ap1: c0=%g", fc, state->a2.ap1.c0);
}

static void cap5_run(struct cap5_state *state, sample_t s, sample_t *lp, sample_t *hp)
{
	sample_t a1 = ap2_run(&state->a1, s);
	sample_t a2 = ap3_run(&state->a2, s);
	*lp = (a1+a2)*0.5;
	*hp = (a1-a2)*0.5;
}

static double fb_bands[5]    = { 250.0, 500.0, 1000.0, 2000.0, 4000.0 };
static int    fb_ap_bands[6] = { 3, 4, 1, 0, 0, 4 };

static void filter_bank_init(struct filter_bank *fb, double fs)
{
	int i;
	for (i = 0; i < 5; ++i)
		cap5_init(&fb->b[i], fs, fb_bands[i]);
	for (i = 0; i < 6; ++i)
		fb->ap[i] = fb->b[fb_ap_bands[i]].a1;
}

static void filter_bank_run(struct filter_bank *fb, sample_t s)
{
	cap5_run(&fb->b[2], s, &fb->s[2], &fb->s[3]);  /* split in the middle (band 2) */
	fb->s[2] = ap2_run(&fb->ap[0], fb->s[2]);      /* band 3 ap */
	fb->s[2] = ap2_run(&fb->ap[1], fb->s[2]);      /* band 4 ap */
	fb->s[3] = ap2_run(&fb->ap[2], fb->s[3]);      /* band 1 ap */
	fb->s[3] = ap2_run(&fb->ap[3], fb->s[3]);      /* band 0 ap */

	cap5_run(&fb->b[1], fb->s[2], &fb->s[1], &fb->s[2]);  /* split at band 1 */
	fb->s[2] = ap2_run(&fb->ap[4], fb->s[2]);             /* band 0 ap */

	cap5_run(&fb->b[3], fb->s[3], &fb->s[3], &fb->s[4]); /* split at band 3 */
	fb->s[3] = ap2_run(&fb->ap[5], fb->s[3]);            /* band 4 ap */

	cap5_run(&fb->b[0], fb->s[1], &fb->s[0], &fb->s[1]); /* split at band 0 */

	cap5_run(&fb->b[4], fb->s[4], &fb->s[4], &fb->s[5]); /* split at band 4 */
}

#define TO_DEGREES(x) ((x)*M_1_PI*180.0)
#define TIME_TO_FREQ(x) (350.0/(x))  /* f0 ~= 0.35/tr */
#define TIME_TO_FRAMES(x, fs) ((x)/1000.0 * (fs))
#define ANGLE(n, d, expr) (((n) == 0 && (d) == 0) ? M_PI_4 : ((d) == 0) ? M_PI_2 : atan(expr))
#define CALC_LR(n, d, expr) (M_PI_4 - ANGLE(n, d, expr))
#define CALC_CS(n, d, expr) (ANGLE(n, d, expr) - M_PI_4)
#define CALC_DIR_BOOST(f, s, n) ((1.0/sqrt(1.0+pow((f)*(s), 2.0))) * (1.0/(n)) - 1.0)

sample_t * matrix4_effect_run(struct effect *e, ssize_t *frames, sample_t *ibuf, sample_t *obuf)
{
	ssize_t i, k, oframes = 0;
	double fl_boost = 0, fr_boost = 0;
	struct matrix4_state *state = (struct matrix4_state *) e->data;

	const double norm_mult = (state->disable) ? 1.0 : state->norm_mult;
	const double surr_mult = (state->disable) ? 0.0 : state->surr_mult;

	for (i = 0; i < *frames; ++i) {
		double f_boost_norm = 0;
		sample_t out_ls = 0, out_rs = 0;

		const double s0 = (ibuf) ? ibuf[i*e->istream.channels + state->c0] : 0.0;
		const double s1 = (ibuf) ? ibuf[i*e->istream.channels + state->c1] : 0.0;

		const double s0_d = state->bufs[state->c0][state->p];
		const double s1_d = state->bufs[state->c1][state->p];

		filter_bank_run(&state->fb[0], s0);
		filter_bank_run(&state->fb[1], s1);

		fl_boost = fr_boost = 0;

		for (k = 0; k < 5; ++k) {
			struct matrix4_band *band = &state->band[k];
			double lr, cs, lr_adapt, cs_adapt;

			const double s0_bp = state->fb[0].s[k+1];
			const double s1_bp = state->fb[1].s[k+1];
			const double s0_d_bp = state->fb_buf[0][k+1][state->p];
			const double s1_d_bp = state->fb_buf[1][k+1][state->p];

			const double lr_ord = CALC_LR(s0_d_bp, s1_d_bp, fabs(s0_d_bp)/fabs(s1_d_bp));
			const double cs_ord = CALC_CS(s0_d_bp+s1_d_bp, s0_d_bp-s1_d_bp, fabs(s0_d_bp+s1_d_bp)/fabs(s0_d_bp-s1_d_bp));

			const double l_pwr = s0_bp*s0_bp;
			const double r_pwr = s1_bp*s1_bp;
			const double sum_pwr = (s0_bp+s1_bp)*(s0_bp+s1_bp);
			const double diff_pwr = (s0_bp-s1_bp)*(s0_bp-s1_bp);

			const double l_pwr_env = ewma_run(&band->r_fast[0], l_pwr);
			const double r_pwr_env = ewma_run(&band->r_fast[1], r_pwr);
			const double sum_pwr_env = ewma_run(&band->r_fast[2], sum_pwr);
			const double diff_pwr_env = ewma_run(&band->r_fast[3], diff_pwr);

			const double l_adapt = l_pwr_env - ewma_run_set_max(&band->adapt[0], l_pwr_env);
			const double r_adapt = r_pwr_env - ewma_run_set_max(&band->adapt[1], r_pwr_env);
			const double sum_adapt = sum_pwr_env - ewma_run_set_max(&band->adapt[2], sum_pwr_env);
			const double diff_adapt = diff_pwr_env - ewma_run_set_max(&band->adapt[3], diff_pwr_env);

			lr_adapt = CALC_LR(l_adapt, r_adapt, sqrt(l_adapt/r_adapt));
			cs_adapt = CALC_CS(sum_adapt, diff_adapt, sqrt(sum_adapt/diff_adapt));

			ewma_run(&band->norm[0], l_adapt);
			ewma_run(&band->norm[1], r_adapt);
			ewma_run(&band->norm[2], l_pwr_env*(1.0/(NORM_MAX)));
			ewma_run(&band->norm[3], r_pwr_env*(1.0/(NORM_MAX)));
			const double l_norm_div_adapt = sqrt(pow(ewma_get_last(&band->norm[0]), 2.0) + pow(ewma_get_last(&band->norm[1])*(NORM_CROSSFEED), 2.0));
			const double r_norm_div_adapt = sqrt(pow(ewma_get_last(&band->norm[1]), 2.0) + pow(ewma_get_last(&band->norm[0])*(NORM_CROSSFEED), 2.0));
			const double l_norm_div_env = sqrt(pow(ewma_get_last(&band->norm[2]), 2.0) + pow(ewma_get_last(&band->norm[3])*(NORM_CROSSFEED), 2.0));
			const double r_norm_div_env = sqrt(pow(ewma_get_last(&band->norm[3]), 2.0) + pow(ewma_get_last(&band->norm[2])*(NORM_CROSSFEED), 2.0));
			const double l_norm_div = MAXIMUM(l_norm_div_adapt, l_norm_div_env);
			const double r_norm_div = MAXIMUM(r_norm_div_adapt, r_norm_div_env);
			const double l_adapt_norm = (l_norm_div > 0.0) ? l_adapt / l_norm_div : (l_adapt <= 0.0) ? 0.0 : 100.0;
			const double r_adapt_norm = (r_norm_div > 0.0) ? r_adapt / r_norm_div : (r_adapt <= 0.0) ? 0.0 : 100.0;
			const double l_adapt_norm_sm = ewma_run(&band->smooth[0], l_adapt_norm);
			const double r_adapt_norm_sm = ewma_run(&band->smooth[1], r_adapt_norm);

			const double l_event = l_adapt_norm_sm - ewma_run(&band->r_slow[0], l_adapt_norm_sm);
			const double r_event = r_adapt_norm_sm - ewma_run(&band->r_slow[1], r_adapt_norm_sm);

			if (band->ev == EVENT_STATE_NONE
					&& (l_event > EVENT_THRESH || r_event > EVENT_THRESH)) {
				band->ev = EVENT_STATE_SAMPLE;
				band->ev_flags = 0;
				band->ev_flags |= (l_event >= r_event) ? EVENT_FLAG_L : 0;
				band->ev_flags |= (r_event >= l_event) ? EVENT_FLAG_R : 0;
				band->t = 0;
				ewma_set(&band->avg[0], lr_adapt);
				ewma_set(&band->avg[1], cs_adapt);
			}

			switch (band->ev) {
			case EVENT_STATE_SAMPLE:
				ewma_run(&band->avg[0], lr_adapt);
				ewma_run(&band->avg[1], cs_adapt);
				lr = ewma_run(&band->drift[0], lr_ord);
				cs = ewma_run(&band->drift[1], cs_ord);
				if (band->t >= state->ev_sample_frames) {
					band->ev = EVENT_STATE_HOLD;
					band->t = 0;
					ewma_set(&band->drift[2], lr);
					ewma_set(&band->drift[3], cs);
					const double abs_lr_ev = fabs(ewma_get_last(&band->avg[0]));
					const double abs_cs_ev = fabs(ewma_get_last(&band->avg[1]));
					if (abs_lr_ev+abs_cs_ev > M_PI_4*1.01) band->ev_flags |= EVENT_FLAG_USE_ORD;
					if (band->ev_flags & EVENT_FLAG_USE_ORD) ++band->ord_count;
					else ++band->diff_count;
					/* LOG_FMT(LL_VERBOSE, "%s: event: type: %4s; lr: %+06.2f°; cs: %+06.2f°",
							e->name, (band->ev_flags & EVENT_FLAG_USE_ORD) ? "ord" : "diff",
							TO_DEGREES(ewma_get_last(&band->avg[0])), TO_DEGREES(ewma_get_last(&band->avg[1]))); */
				}
				else ++band->t;
				break;
			case EVENT_STATE_HOLD:
				if (band->ev_flags & EVENT_FLAG_USE_ORD) {
					lr_adapt = lr = ewma_run(&band->drift[0], lr_ord);
					cs_adapt = cs = ewma_run(&band->drift[1], cs_ord);
				}
				else {
					lr_adapt = lr = ewma_set(&band->drift[0], ewma_run(&band->drift[2], ewma_get_last(&band->avg[0])));
					cs_adapt = cs = ewma_set(&band->drift[1], ewma_run(&band->drift[3], ewma_get_last(&band->avg[1])));
				}
				if ((band->ev_flags & EVENT_FLAG_L && l_adapt_norm_sm <= EVENT_END_THRESH)
						|| (band->ev_flags & EVENT_FLAG_R && r_adapt_norm_sm <= EVENT_END_THRESH)) {
					band->ev_flags |= EVENT_FLAG_END;
				}
				if ((band->t >= state->ev_min_hold_frames && band->ev_flags & EVENT_FLAG_END)
						|| band->t >= state->ev_max_hold_frames) {
					if (band->t < state->ev_max_hold_frames) ++band->early_count;
					band->ev = EVENT_STATE_NONE;
					band->t = 0;
				}
				else ++band->t;
				break;
			default:
				lr = ewma_run(&band->drift[0], lr_ord);
				cs = ewma_run(&band->drift[1], cs_ord);
				lr_adapt = cs_adapt = 0;
				++band->t;
			}

			double abs_lr = fabs(lr);
			double abs_cs = fabs(cs);
			if (abs_lr+abs_cs > M_PI_4) {
				const double norm = M_PI_4 / (abs_lr+abs_cs);
				lr *= norm;
				cs *= norm;
				abs_lr *= norm;
				abs_cs *= norm;
			}

			sample_t lsl_m, lsr_m, rsl_m, rsr_m;

			if (state->do_dir_boost) {
				if (cs >= 0.0) {
					const double dir_boost = CALC_DIR_BOOST(cos(2.0*(abs_lr+cs)), surr_mult, norm_mult);
					band->fl_boost = (lr > 0.0) ? dir_boost*(1.0-sin(2.0*lr)) : dir_boost;
					band->fr_boost = (lr < 0.0) ? dir_boost*(1.0-sin(-2.0*lr)) : dir_boost;
				}
				else if (cs >= -M_PI_4/2) {
					const double dir_boost = CALC_DIR_BOOST(1.0-(1.0-cos(2.0*lr))*cos(4.0*cs), surr_mult, norm_mult);
					band->fl_boost = (lr > 0.0) ? dir_boost*(1.0-sin(2.0*lr)) : dir_boost;
					band->fr_boost = (lr < 0.0) ? dir_boost*(1.0-sin(-2.0*lr)) : dir_boost;
				}
				else {
					band->fl_boost = 0.0;
					band->fr_boost = 0.0;
				}
			}
			else {
				band->fl_boost = 0.0;
				band->fr_boost = 0.0;
			}
			const double d_sum_pwr_env = ewma_run(&state->band_d_sum_pwr_env[k+1], (s0_d_bp+s1_d_bp)*(s0_d_bp+s1_d_bp));
			fl_boost += band->fl_boost * band->fl_boost * d_sum_pwr_env;
			fr_boost += band->fr_boost * band->fr_boost * d_sum_pwr_env;
			f_boost_norm += d_sum_pwr_env;

			/* The matrix coefficients for the surround channels are from
			   "Multichannel matrix surround decoders for two-eared listeners" by
			   David Griesinger (http://www.davidgriesinger.com/sur.pdf). I've
			   corrected gsl so there is full cancellation when |lr|+|cs|=45°.
			*/
			const double gl = (cos(M_PI/4-abs_lr)-sin(M_PI/4-abs_lr))/cos(M_PI/4-abs_lr);
			const double gsl = gl*gl;
			if (cs >= 0.0) {
				const double cf = cos(cs)+sin(cs);
				const double gc = 2.0*sin(cs)/(cos(cs)+sin(cs));
				if (lr <= 0.0) {
					lsl_m = cf*(1.0-gsl-0.5*gc);
					lsr_m = cf*(-0.5*gc-gl);
					rsl_m = -sin(cs);
					rsr_m = cos(cs);
				}
				else {
					lsl_m = cos(cs);
					lsr_m = -sin(cs);
					rsl_m = cf*(-0.5*gc-gl);
					rsr_m = cf*(1.0-gsl-0.5*gc);
				}
			}
			else {
				if (lr <= 0.0) {
					if (cs > -M_PI_4/2) {
						lsl_m = 1.0-gsl*(1.0+sin(4.0*cs));
						lsr_m = -gl*cos(4.0*cs);
					}
					else {
						lsl_m = 1.0;
						lsr_m = 0.0;
					}
					rsl_m = 0.0;
					rsr_m = 1.0;
				}
				else {
					lsl_m = 1.0;
					lsr_m = 0.0;
					if (cs > -M_PI_4/2) {
						rsl_m = -gl*cos(4.0*cs);
						rsr_m = 1.0-gsl*(1.0+sin(4.0*cs));
					}
					else {
						rsl_m = 0.0;
						rsr_m = 1.0;
					}
				}
			}

			const sample_t ls_m_pwr = sqrt(lsl_m*lsl_m + lsr_m*lsr_m);
			const sample_t rs_m_pwr = sqrt(rsl_m*rsl_m + rsr_m*rsr_m);
			lsl_m /= ls_m_pwr;
			lsr_m /= ls_m_pwr;
			rsl_m /= rs_m_pwr;
			rsr_m /= rs_m_pwr;

			out_ls += (s0_d_bp*lsl_m + s1_d_bp*lsr_m) * norm_mult * surr_mult;
			out_rs += (s0_d_bp*rsl_m + s1_d_bp*rsr_m) * norm_mult * surr_mult;

			if (k == 0) {
				const double eb_s0_d_bp = state->fb_buf[0][0][state->p];
				const double eb_s1_d_bp = state->fb_buf[1][0][state->p];
				const double eb_sum_d_bp = biquad(&state->band0_hp, eb_s0_d_bp+eb_s1_d_bp);
				out_ls += (eb_s0_d_bp*lsl_m + eb_s1_d_bp*lsr_m) * norm_mult * surr_mult;
				out_rs += (eb_s0_d_bp*rsl_m + eb_s1_d_bp*rsr_m) * norm_mult * surr_mult;
				const double eb_d_sum_pwr_env = ewma_run(&state->band_d_sum_pwr_env[0], eb_sum_d_bp*eb_sum_d_bp);
				fl_boost += band->fl_boost * band->fl_boost * eb_d_sum_pwr_env;
				fr_boost += band->fr_boost * band->fr_boost * eb_d_sum_pwr_env;
				f_boost_norm += eb_d_sum_pwr_env;
			}

			band->lr = lr;
			band->cs = cs;
			band->lr_adapt = lr_adapt;
			band->cs_adapt = cs_adapt;
		}

		/* FIXME: Weighted RMS average based on power envelope... maybe there's a better way? */
		if (f_boost_norm > 0.0) {
			fl_boost = sqrt(fl_boost / f_boost_norm);
			fr_boost = sqrt(fr_boost / f_boost_norm);
		}
		else {
			fl_boost = 0.0;
			fr_boost = 0.0;
		}

		sample_t ll_m, lr_m, rl_m, rr_m;
		ll_m = 1.0 + fl_boost;
		lr_m = 0.0;
		rl_m = 0.0;
		rr_m = 1.0 + fr_boost;

		const sample_t out_l = (s0_d*ll_m + s1_d*lr_m) * norm_mult;
		const sample_t out_r = (s0_d*rl_m + s1_d*rr_m) * norm_mult;

		for (k = 0; k < 6; ++k) {
			state->fb_buf[0][k][state->p] = state->fb[0].s[k];
			state->fb_buf[1][k][state->p] = state->fb[1].s[k];
		}

		if (state->has_output) {
			for (k = 0; k < e->istream.channels; ++k) {
				if (k == state->c0)
					obuf[oframes*e->ostream.channels + k] = out_l;
				else if (k == state->c1)
					obuf[oframes*e->ostream.channels + k] = out_r;
				else
					obuf[oframes*e->ostream.channels + k] = state->bufs[k][state->p];
				state->bufs[k][state->p] = (ibuf) ? ibuf[i*e->istream.channels + k] : 0.0;
			}
			obuf[oframes*e->ostream.channels + k + 0] = out_ls;
			obuf[oframes*e->ostream.channels + k + 1] = out_rs;
			++oframes;
		}
		else {
			for (k = 0; k < e->istream.channels; ++k) {
				#ifdef SYMMETRIC_IO
					obuf[oframes*e->ostream.channels + k] = 0.0;
				#endif
				state->bufs[k][state->p] = (ibuf) ? ibuf[i*e->istream.channels + k] : 0.0;
			}
			#ifdef SYMMETRIC_IO
				obuf[oframes*e->ostream.channels + k + 0] = 0.0;
				obuf[oframes*e->ostream.channels + k + 1] = 0.0;
				++oframes;
			#endif
		}
		state->p = (state->p + 1 >= state->len) ? 0 : state->p + 1;
		if (state->p == 0)
			state->has_output = 1;
	}
	#ifndef LADSPA_FRONTEND
		/* TODO: Implement a proper way for effects to show status lines. */
		if (state->show_status) {
			for (i = 0; i < 5; ++i) {
				fprintf(stderr, "\n%s%s: band %zd: lr: %+06.2f (%+06.2f); cs: %+06.2f (%+06.2f); dir_boost: l:%05.2f r:%05.2f; ord: %zd; diff: %zd; early: %zd\033[K\r",
					e->name, (state->disable) ? " [off]" : "", i,
					TO_DEGREES(state->band[i].lr), TO_DEGREES(state->band[i].lr_adapt), TO_DEGREES(state->band[i].cs), TO_DEGREES(state->band[i].cs_adapt),
					state->band[i].fl_boost, state->band[i].fr_boost, state->band[i].ord_count, state->band[i].diff_count, state->band[i].early_count);
			}
			fprintf(stderr, "\n%s%s: weighted RMS dir_boost: l:%05.2f r:%05.2f\033[K\r",
				e->name, (state->disable) ? " [off]" : "", fl_boost, fr_boost);
			fprintf(stderr, "\033[%zdA", i+1);
		}
	#endif

	*frames = oframes;
	return obuf;
}

ssize_t matrix4_effect_delay(struct effect *e)
{
	struct matrix4_state *state = (struct matrix4_state *) e->data;
	return (state->has_output) ? state->len : state->p;
}

void matrix4_effect_reset(struct effect *e)
{
	int i;
	struct matrix4_state *state = (struct matrix4_state *) e->data;
	state->p = 0;
	state->has_output = 0;
	for (i = 0; i < e->istream.channels; ++i)
		memset(state->bufs[i], 0, state->len * sizeof(sample_t));
	for (i = 0; i < 6; ++i) {
		memset(state->fb_buf[0][i], 0, state->len * sizeof(sample_t));
		memset(state->fb_buf[1][i], 0, state->len * sizeof(sample_t));
	}
}

void matrix4_effect_signal(struct effect *e)
{
	struct matrix4_state *state = (struct matrix4_state *) e->data;
	state->disable = !state->disable;
	LOG_FMT(LL_NORMAL, "%s: %s", e->name, (state->disable) ? "disabled" : "enabled");
}

void matrix4_effect_drain(struct effect *e, ssize_t *frames, sample_t *obuf)
{
	struct matrix4_state *state = (struct matrix4_state *) e->data;
	if (!state->has_output && state->p == 0)
		*frames = -1;
	else {
		if (!state->is_draining) {
			state->drain_frames = state->len;
			state->is_draining = 1;
		}
		if (state->drain_frames > 0) {
			*frames = MINIMUM(*frames, state->drain_frames);
			state->drain_frames -= *frames;
			e->run(e, frames, NULL, obuf);
		}
		else
			*frames = -1;
	}
}

void matrix4_effect_destroy(struct effect *e)
{
	int i;
	struct matrix4_state *state = (struct matrix4_state *) e->data;
	for (i = 0; i < e->istream.channels; ++i)
		free(state->bufs[i]);
	free(state->bufs);
	#ifndef LADSPA_FRONTEND
		if (state->show_status) {
			for (i = 0; i < 6; ++i) fprintf(stderr, "\033[K\n");
			fprintf(stderr, "\033[K\r\033[%dA", i);
		}
	#endif
	free(state);
}

struct effect * matrix4_effect_init(struct effect_info *ei, struct stream_info *istream, char *channel_selector, const char *dir, int argc, char **argv)
{
	struct effect *e;
	struct matrix4_state *state;
	int i, k, n_channels = 0;
	double surr_mult = 0.5;
	char *opt = NULL, *next_opt, *endptr;

	if (argc > 3) {
		LOG_FMT(LL_ERROR, "%s: usage: %s", argv[0], ei->usage);
		return NULL;
	}
	if (argc > 2)
		opt = argv[1];
	if (argc > 1) {
		surr_mult = pow(10.0, strtod(argv[argc-1], &endptr) / 20.0);
		CHECK_ENDPTR(argv[argc-1], endptr, "surround_level", return NULL);
		if (surr_mult > 1.0)
			LOG_FMT(LL_ERROR, "%s: warning: surround_level probably shouldn't be greater than 0", argv[0]);
	}

	if (istream->fs < 32000) {
		LOG_FMT(LL_ERROR, "%s: error: sample rate out of range", argv[0]);
		return NULL;
	}
	for (i = 0; i < istream->channels; ++i)
		if (GET_BIT(channel_selector, i))
			++n_channels;
	if (n_channels != 2) {
		LOG_FMT(LL_ERROR, "%s: error: number of input channels must be 2", argv[0]);
		return NULL;
	}

	e = calloc(1, sizeof(struct effect));
	e->name = ei->name;
	e->istream.fs = e->ostream.fs = istream->fs;
	e->istream.channels = istream->channels;
	e->ostream.channels = istream->channels + 2;
	e->run = matrix4_effect_run;
	e->delay = matrix4_effect_delay;
	e->reset = matrix4_effect_reset;
	e->drain = matrix4_effect_drain;
	e->destroy = matrix4_effect_destroy;
	state = calloc(1, sizeof(struct matrix4_state));

	state->c0 = state->c1 = -1;
	for (i = 0; i < istream->channels; ++i) {
		if (GET_BIT(channel_selector, i)) {
			if (state->c0 == -1)
				state->c0 = i;
			else
				state->c1 = i;
		}
	}
	state->do_dir_boost = 1;
	if (opt) {
		while (*opt != '\0') {
			next_opt = isolate(opt, ',');
			if (*opt == '\0') /* do nothing */;
			else if (strcmp(opt, "show_status")  == 0) state->show_status = 1;
			else if (strcmp(opt, "no_dir_boost") == 0) state->do_dir_boost = 0;
			else if (strcmp(opt, "signal")       == 0) e->signal = matrix4_effect_signal;
			else {
				LOG_FMT(LL_ERROR, "%s: error: unrecognized option: %s", argv[0], opt);
				goto opt_fail;
			}
			opt = next_opt;
		}
	}
	filter_bank_init(&state->fb[0], istream->fs);
	filter_bank_init(&state->fb[1], istream->fs);
	for (k = 0; k < 5; ++k) {
		for (i = 0; i < 4; ++i) ewma_init(&state->band[k].r_fast[i], istream->fs, TIME_TO_FREQ(RISE_TIME_FAST));
		for (i = 0; i < 2; ++i) ewma_init(&state->band[k].r_slow[i], istream->fs, TIME_TO_FREQ(RISE_TIME_SLOW));
		for (i = 0; i < 4; ++i) ewma_init(&state->band[k].adapt[i], istream->fs, TIME_TO_FREQ(ADAPT_TIME));
		for (i = 0; i < 4; ++i) ewma_init(&state->band[k].norm[i], istream->fs, TIME_TO_FREQ(NORM_TIME));
		for (i = 0; i < 2; ++i) ewma_init(&state->band[k].smooth[i], istream->fs, TIME_TO_FREQ(RISE_TIME_FAST));
		for (i = 0; i < 2; ++i) ewma_init(&state->band[k].avg[i], istream->fs, TIME_TO_FREQ(EVENT_SAMPLE_TIME));
		for (i = 0; i < 2; ++i) ewma_init(&state->band[k].drift[i], istream->fs, TIME_TO_FREQ(ADAPT_TIME));
		for (i = 2; i < 4; ++i) ewma_init(&state->band[k].drift[i], istream->fs, TIME_TO_FREQ(RISE_TIME_FAST));
	}
	for (k = 0; k < 6; ++k)
		ewma_init(&state->band_d_sum_pwr_env[k], istream->fs, TIME_TO_FREQ(RISE_TIME_FAST));
	biquad_init_using_type(&state->band0_hp, BIQUAD_HIGHPASS, istream->fs, 38.0, 0.5, 0, 0, BIQUAD_WIDTH_Q);

	state->len = lround(TIME_TO_FRAMES(RISE_TIME_FAST + EVENT_SAMPLE_TIME, istream->fs));
	state->bufs = calloc(istream->channels, sizeof(sample_t *));
	for (i = 0; i < istream->channels; ++i)
		state->bufs[i] = calloc(state->len, sizeof(sample_t));
	for (i = 0; i < 6; ++i) {
		state->fb_buf[0][i] = calloc(state->len, sizeof(sample_t));
		state->fb_buf[1][i] = calloc(state->len, sizeof(sample_t));
	}
	state->surr_mult = surr_mult;
	state->norm_mult = 1.0 / sqrt(1.0 + surr_mult*surr_mult);
	state->ev_sample_frames = TIME_TO_FRAMES(EVENT_SAMPLE_TIME, istream->fs);
	state->ev_max_hold_frames = TIME_TO_FRAMES(EVENT_MAX_HOLD_TIME, istream->fs);
	state->ev_min_hold_frames = TIME_TO_FRAMES(EVENT_MIN_HOLD_TIME, istream->fs);
	e->data = state;
	return e;

	opt_fail:
	free(state);
	free(e);
	return NULL;
}
