#define _ncCheck(X) ncerr = X ; \
if (ncerr /= nf90_noerr) then ; \
print *, trim(nf90_strerror(ncerr)) ; \
error = -1 ; \
return ; \
endif
