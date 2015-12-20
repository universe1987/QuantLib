======================================================
General Todos
======================================================

We have errors from boost that an unqualified sqrt call (e.g. L41 in error_of_mean) is ambiguous. Analyze this. In tmp/cppad_error_of_mean.cpp is an example, where the error is *not* reproduced though. Why is it happening in QL then? Fix: In boost (in folder Peter only), revert this !!! Also in QL I replaced boost::normal in KahaleSmileSection with QL's implementation for the same reason.

======================================================
Update qrm_cvarisk_adjoint on 15-Dec-2015 from initial
fork point to current master (peter), which is similar
to the official QuantLib release 1.7
======================================================

The initial fork point was identified as 8ac8fe2b (quantlib-old, i.e. before the modules were separated) by doing

git merge-base --all master adjoint

The branch qrm_cvarisk_adjoint was then updated from this point to the current master (peter) 0416fed7dd (quantlib-old, 706f3c6b473383 after separation of modules), release 1.7 was just released in the official QuantLib project.

The status in qrm_cvarisk_adjoint before the update is c07e932b (quantlib-old, 4e408d0e06f48ea after separation of modules)

Notes:
- did not upgrade Examples
- did not upgrade msvc projects, dev (not supported anyway)
- did not insert certain proprietary changes from peter's master, in particular
-- last period day counter for fixed rate coupons
-- new inflation coupon types
-- monte carlo multithreading contribution
-- beta eta model implementation
-- quasi gaussian model implementation draft
-- matrix access extensions (a la boost)
-- Libor market model extensions
-- Timegrid extensions

Todos:
- templatize again (since overwritten with updates --- it seemed less work to convert again than to manually merge the changes):
-- gaussian1dmodel, gsr, gsrprocess
-- experimental/volatility -- some classes needed to be excluded from build (convert them later and add them again)

Lessons learned for the next upgrade:
- document the changes in couponpricers, coupons (pricer are moved, sub-files are created ...)
- build tools that replace forward declarations like class YieldTermStructure; by template<class T> YieldTermStructure_t; typedef YieldTermStructure_t<Real> YieldTermStructure;

