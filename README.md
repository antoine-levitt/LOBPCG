# LOBPCG

Implementation of the LOBPCG algorithm, seeking optimal performance without compromising on stability. The code is written for readability, and is not fully optimized yet (ie it allocates a lot of memory). The most time-consuming parts (BLAS3 and matvecs) should be OK though. It follows the scheme of Hetmaniuk & Lehoucq (see also the refinements in Duersch et. al.), with the following modifications:

* Cholesky factorizations are used instead of the eigenvalue decomposition of Stathopolous & Wu for the orthogonalization. Cholesky orthogonalization has an unwarranted bad reputation: when applied twice to a matrix X with κ(X) <~ sqrt(1/ε), where ε is the machine epsilon, it will produce a set of very orthogonal vectors, just as the eigendecomposition-based method. It can fail when κ >~ sqrt(1/ε), but that can be fixed by shifting the overlap matrix. This is 100% reliable while being much faster than eigendecompositions. 
* Termination criteria for the orthogonalization are based on cheap estimates instead of costly explicit checks.
* The default tolerances are very tight, yielding a very stable algorithm. This can be tweaked (see keyword `ortho_tol`) to compromise on stability for less Cholesky factorizations.
* Implicit product updates (reuse of matvecs) are performed whenever it is safe to do so, ie only with orthogonal transformations. An exception is the B matrix reuse, which seems to be OK even with very badly conditioned B matrices. The code is easy to modify if this is not the case.
* The locking is performed carefully to ensure minimal impact on the other eigenvectors (which is not the case if performed naively)
