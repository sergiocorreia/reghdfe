
(* TO DO *)


Note, I also need to go into detail about a few things:

1. Why is there no constant? Because it made the code messy everywhere, and because in practice it's just absorbed by the FEs (and you can recover it from them, by taking the average of the predict,d)

2. There are four R2s. Overall/within, and standard/adjusted. The formula for adjusted/within is different from the one in xtreg/xtivreg2 which is IMHO wrong. See forum post notes, see montecarlo sim, etc.