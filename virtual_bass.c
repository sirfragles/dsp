#include <stdlib.h>
#include <math.h>
#include "virtual_bass.h"
#include "biquad.h"
#include "util.h"

struct virtual_bass_state {
	sample_t level, envelope, time_const;
	struct biquad_state lp[3];
	struct biquad_state hp[3];
};

static inline sample_t harmonic_func(sample_t s, sample_t e, sample_t l)
{
	e = (e > 0.0) ? MINIMUM(1.0/e, 1000.0) : 1000.0;
	sample_t e2 = e*e; sample_t e3 = e2*e;
	sample_t s2 = s*s; sample_t s3 = s2*s; sample_t s4 = s3*s;
	return e*s2*l + e2*s3*l*0.5 + e3*s4*l*0.25;
}

sample_t * virtual_bass_effect_run(struct effect *e, ssize_t *frames, sample_t *ibuf, sample_t *obuf)
{
	sample_t s, lf, hf;
	ssize_t samples = *frames * e->ostream.channels, i, k;
	struct virtual_bass_state **state = (struct virtual_bass_state **) e->data;
	for (i = 0; i < samples; i += e->ostream.channels) {
		for (k = 0; k < e->ostream.channels; ++k) {
			if (state[k]) {
				s = ibuf[i + k];
				lf = biquad(&state[k]->lp[1], biquad(&state[k]->lp[0], s));
				hf = biquad(&state[k]->hp[1], biquad(&state[k]->hp[0], s));
				state[k]->envelope = MAXIMUM(state[k]->envelope, fabs(lf));
				ibuf[i + k] = hf + lf + biquad(&state[k]->hp[2], biquad(&state[k]->lp[2], harmonic_func(lf, state[k]->envelope, state[k]->level)));
				state[k]->envelope *= state[k]->time_const;
			}
		}
	}
	return ibuf;
}

void virtual_bass_effect_reset(struct effect *e)
{
	int i;
	struct virtual_bass_state **state = (struct virtual_bass_state **) e->data;
	for (i = 0; i < e->ostream.channels; ++i) {
		if (state[i]) {
			biquad_reset(&state[i]->lp[0]);
			biquad_reset(&state[i]->lp[1]);
			biquad_reset(&state[i]->lp[2]);
			biquad_reset(&state[i]->hp[0]);
			biquad_reset(&state[i]->hp[1]);
			biquad_reset(&state[i]->hp[2]);
		}
	}
}

void virtual_bass_effect_destroy(struct effect *e)
{
	int i;
	struct virtual_bass_state **state = (struct virtual_bass_state **) e->data;
	for (i = 0; i < e->ostream.channels; ++i)
		free(state[i]);
	free(state);
}

struct effect * virtual_bass_effect_init(struct effect_info *ei, struct stream_info *istream, char *channel_selector, const char *dir, int argc, char **argv)
{
	int i;
	double cutoff, level;
	struct virtual_bass_state **state;
	struct effect *e;
	char *endptr;

	if (argc < 2 || argc > 3) {
		LOG_FMT(LL_ERROR, "%s: usage: %s", argv[0], ei->usage);
		return NULL;
	}

	cutoff = parse_freq(argv[1], &endptr);
	CHECK_ENDPTR(argv[1], endptr, "cutoff", return NULL);
	CHECK_FREQ(cutoff, istream->fs/5.0, "cutoff", return NULL);
	if (argc == 3) {
		level = strtod(argv[2], &endptr);
		CHECK_ENDPTR(argv[2], endptr, "level", return NULL);
		CHECK_RANGE(level >= 0.0, "level", return NULL);
	}
	else
		level = 1.0;

	e = calloc(1, sizeof(struct effect));
	e->name = ei->name;
	e->istream.fs = e->ostream.fs = istream->fs;
	e->istream.channels = e->ostream.channels = istream->channels;
	e->run = virtual_bass_effect_run;
	e->reset = virtual_bass_effect_reset;
	e->destroy = virtual_bass_effect_destroy;
	state = calloc(istream->channels, sizeof(struct virtual_bass_state *));
	for (i = 0; i < istream->channels; ++i) {
		if (GET_BIT(channel_selector, i)) {
			state[i] = calloc(1, sizeof(struct virtual_bass_state));
			state[i]->level = level;
			state[i]->time_const = pow(0.5, cutoff / 100.0 / istream->fs);
			biquad_init_using_type(&state[i]->lp[0], BIQUAD_LOWPASS, istream->fs, cutoff, 0.5, 0.0, 0.0, BIQUAD_WIDTH_Q);
			biquad_init_using_type(&state[i]->lp[1], BIQUAD_LOWPASS, istream->fs, cutoff, 1.0, 0.0, 0.0, BIQUAD_WIDTH_Q);
			biquad_init_using_type(&state[i]->lp[2], BIQUAD_LOWPASS, istream->fs, cutoff * 2.5, 0.707106781187, 0.0, 0.0, BIQUAD_WIDTH_Q);
			biquad_init_using_type(&state[i]->hp[0], BIQUAD_HIGHPASS, istream->fs, cutoff, 0.5, 0.0, 0.0, BIQUAD_WIDTH_Q);
			biquad_init_using_type(&state[i]->hp[1], BIQUAD_HIGHPASS, istream->fs, cutoff, 1.0, 0.0, 0.0, BIQUAD_WIDTH_Q);
			biquad_init_using_type(&state[i]->hp[2], BIQUAD_HIGHPASS, istream->fs, cutoff * 0.75, 0.707106781187, 0.0, 0.0, BIQUAD_WIDTH_Q);
		}
	}
	e->data = state;
	return e;
}
