(UPDATE THIS)

# Quick wishlist for possible future additions:

* Expand the hdfe.ado algorithm, so that we can switch between the current algorithm, the one in Somaini & Wolak (for the low-dimensional FEs), and the one by Mittag (c.g., for an initial acceleration).
* If what we want is *speed*, nothing would gain rewriting the core algorithm in C and even CUDA/GPU. However, that is extremely low on the list of things I'm interested in doing (and I also don't know CUDA).
* Add core() support for OSX/Linux. One problem is that there is no easy way to run winexec on terminal. Maybe we can switch to -parallel- entirely, but that may require some rewriting.
* Maybe we can implement a test for the FEs under robust/cluster with the robustified hausman. See explanation of xtoverid (Hansen J)
* Add a more thorough discussion on the possible identification issues. Maybe as a PDF
* Implement -bootstrap- option in DoF estimation
* Calculate exact DoF adjustment for 3+ HDFEs (note: not a problem with cluster VCE when one FE is nested within the cluster). Note: Paulo Guimaraes has worked on this, as well as Simen Gaure with *lfe* (and a few emails received have discussions on the topic).
* Allow a simpler way to replicate -suest-
