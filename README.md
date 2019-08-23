# DanteJulia
Finite volume MHD simulation with structured mesh. This is rewritten from the Matlab version, with improved performance and capabilities.

One issue I encountered is using @view for face values. This greatly slows down the flux calculations because of heap allocated memories. As suggested by Roger on Julia's Chinese forum, one workaround is to use the unsafeArrays.jl package to allocate @view memory on the stack. This is worth trying because for the current implementation, LState_XV and RState_XV are just shifts of the original array State_GV. Ideally there is no need to copy the data: using pointers/views should be enough.

Let me create a simple scenario to deal with the problem and find out a solution. The package UNsafeArray.jl is worth trying.