#ifndef _DOLBY_ATMOS_512_H
#define _DOLBY_ATMOS_512_H

#include "dsp.h"
#include "effect.h"

struct effect * dolby_atmos_512_effect_init(struct effect_info *, struct stream_info *, char *, const char *, int, char **);

#endif
