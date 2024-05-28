#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gain.h"
#include "util.h"

#define MAX_OPERATIONS 10

typedef enum {
    OPERATION_GAIN,
    OPERATION_MULT,
    OPERATION_ADD,
    OPERATION_PHASE_REVERSE
} operation_type;

struct operation {
    operation_type type;
    sample_t value;
};

struct gain_state {
    int channel;
    int num_operations;
    struct operation operations[MAX_OPERATIONS];
};

sample_t * gain_effect_run(struct effect *e, ssize_t *frames, sample_t *ibuf, sample_t *obuf)
{
    ssize_t i, k, samples = *frames * e->ostream.channels;
    struct gain_state *state = (struct gain_state *) e->data;

    for (int op_index = 0; op_index < state->num_operations; ++op_index) {
        struct operation *op = &state->operations[op_index];
        if (state->channel == -1) {
            for (i = 0; i < samples; i += e->ostream.channels) {
                for (k = 0; k < e->ostream.channels; ++k) {
                    if (GET_BIT(e->channel_selector, k)) {
                        switch (op->type) {
                            case OPERATION_GAIN:
                            case OPERATION_MULT:
                            case OPERATION_PHASE_REVERSE:
                                ibuf[i + k] *= op->value;
                                break;
                            case OPERATION_ADD:
                                ibuf[i + k] += op->value;
                                break;
                        }
                    }
                }
            }
        } else {
            for (i = state->channel; i < samples; i += e->ostream.channels) {
                switch (op->type) {
                    case OPERATION_GAIN:
                    case OPERATION_MULT:
                    case OPERATION_PHASE_REVERSE:
                        ibuf[i] *= op->value;
                        break;
                    case OPERATION_ADD:
                        ibuf[i] += op->value;
                        break;
                }
            }
        }
    }
    return ibuf;
}

void gain_effect_plot(struct effect *e, int i)
{
    struct gain_state *state = (struct gain_state *) e->data;
    int k;
    for (int op_index = 0; op_index < state->num_operations; ++op_index) {
        struct operation *op = &state->operations[op_index];
        if (state->channel == -1) {
            for (k = 0; k < e->ostream.channels; ++k) {
                if (GET_BIT(e->channel_selector, k))
                    printf("H%d_%d(f)=%.15e\n", k, i, 20 * log10(fabs(op->value)));
                else
                    printf("H%d_%d(f)=0\n", k, i);
            }
        } else {
            for (k = 0; k < e->ostream.channels; ++k) {
                if (k == state->channel)
                    printf("H%d_%d(f)=%.15e\n", k, i, 20 * log10(fabs(op->value)));
                else
                    printf("H%d_%d(f)=0\n", k, i);
            }
        }
    }
}

void gain_effect_destroy(struct effect *e)
{
    free(e->data);
    free(e->channel_selector);
}

struct effect * gain_effect_init(struct effect_info *ei, struct stream_info *istream, char *channel_selector, const char *dir, int argc, char **argv)
{
    struct effect *e;
    struct gain_state *state;
    int channel = -1;
    char *endptr;

    if (argc < 2 || argc > (MAX_OPERATIONS + 1)) {
        LOG_FMT(LL_ERROR, "%s: usage: %s", argv[0], ei->usage);
        return NULL;
    }
    if (argc >= 3) {
        channel = strtol(argv[1], &endptr, 10);
        CHECK_ENDPTR(argv[1], endptr, "channel", return NULL);
        CHECK_RANGE(channel >= 0 && channel < istream->channels, "channel", return NULL);
    }

    e = calloc(1, sizeof(struct effect));
    e->name = ei->name;
    e->istream.fs = e->ostream.fs = istream->fs;
    e->istream.channels = e->ostream.channels = istream->channels;
    e->channel_selector = NEW_SELECTOR(istream->channels);
    COPY_SELECTOR(e->channel_selector, channel_selector, istream->channels);

    state = calloc(1, sizeof(struct gain_state));
    state->channel = channel;
    state->num_operations = 0;

    for (int arg_index = (channel == -1) ? 1 : 2; arg_index < argc; ++arg_index) {
        if (strcmp(argv[arg_index], "phase_revers") == 0) {
            state->operations[state->num_operations].type = OPERATION_PHASE_REVERSE;
            state->operations[state->num_operations].value = -1.0;
        } else {
            double value = strtod(argv[arg_index], &endptr);
            CHECK_ENDPTR(argv[arg_index], endptr, "value", return NULL);

            if (state->num_operations == 0 && ei->effect_number == GAIN_EFFECT_NUMBER_GAIN) {
                value = pow(10.0, value / 20.0);
                state->operations[state->num_operations].type = OPERATION_GAIN;
            } else if (state->num_operations == 0 && ei->effect_number == GAIN_EFFECT_NUMBER_MULT) {
                state->operations[state->num_operations].type = OPERATION_MULT;
            } else if (state->num_operations == 0 && ei->effect_number == GAIN_EFFECT_NUMBER_ADD) {
                state->operations[state->num_operations].type = OPERATION_ADD;
            } else {
                LOG_FMT(LL_ERROR, "%s: invalid argument: %s", argv[0], argv[arg_index]);
                return NULL;
            }
            state->operations[state->num_operations].value = value;
        }
        state->num_operations++;
    }

    e->run = gain_effect_run;
    e->plot = gain_effect_plot;
    e->destroy = gain_effect_destroy;
    e->data = state;

    return e;
}
