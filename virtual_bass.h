#ifndef _VIRTUAL_BASS_H
#define _VIRTUAL_BASS_H

#include "dsp.h"
#include "effect.h"

struct effect * virtual_bass_effect_init(struct effect_info *, struct stream_info *, char *, const char *, int, char **);

#endif
