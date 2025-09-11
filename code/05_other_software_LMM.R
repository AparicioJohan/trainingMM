

mod_asr <- asreml(fixed = yield ~ block + gen, data = data,  family = asr_gaussian(dispersion = 1))
