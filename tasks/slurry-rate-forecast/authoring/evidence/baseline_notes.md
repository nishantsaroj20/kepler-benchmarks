# Baselines and ablations

Tolerances are an absolute rate error of at most the larger of 0.090 times the latent rate and 0.00015 mol m^{-3} s^{-1}, an absolute consumption error of at most the larger of 0.060 times the latent consumption and 0.20 mol m^{-3}, and an absolute fraction error of at most 0.030. Counts are in `calibration_results.json`.

Constant 0, constant 1, the mean of the latent values, the median of the latent values, and copying bath temperature, gauge pressure, catalyst mass, particle diameter, or stirrer speed into every numeric field are outside the band on all 1539 charges. The container run of a program that writes zeros reported the same 1539 misses.

Leaving out internal diffusion, while keeping the reference parameters and integrating the algebraic two-film rate, leaves 1447 charges outside. Maximum relative rate error 0.9783627348125193. Replacing consumption with the initial rate times the hold, using the reference initial state, leaves 1468 charges outside. Refitting on the spreadsheet rate column misses the rate band on 985 charges. Maximum relative rate error 0.3877605151595768.

`naive_forecast.py` is the reference particle balance run one charge at a time. On two charges it matches the batched reference. On the full map the verifier killed it at 20 seconds, discarded the partial output, and recorded reward 0.
