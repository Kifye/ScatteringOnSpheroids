module elemfunc
    use regime
contains
    function qfloat(x)
        integer x
        real(knd) qfloat
        qfloat = real(x, 16)
        return
    end function qfloat

    function qlog10(x)
        real(knd) x, qlog10
        qlog10 = log10(x)
        return
    end function qlog10

    function cqabs(x)
        complex(knd) x
        real(knd) cqabs
        cqabs = abs(x)
        return
    end function cqabs

    function qabs(x)
        real(knd) x, qabs
        qabs = abs(x)
        return
    end function qabs

    function qsin(x)
        real(knd) x, qsin
        qsin = sin(x)
        return
    end function qsin

    function cqsin(x)
        complex(knd) x, cqsin
        cqsin = sin(x)
        return
    end function cqsin

    function qcos(x)
        real(knd) x, qcos
        qcos = cos(x)
        return
    end function qcos

    function qlog(x)
        real(knd) x, qlog
        qlog = log(x)
        return
    end function qlog

    function cqcos(x)
        complex(knd) x, cqcos
        cqcos = cos(x)
        return
    end function cqcos

    function qsqrt(x)
        real(knd) x, qsqrt
        qsqrt = sqrt(x)
        return
    end function qsqrt

    function cqsqrt(x)
        complex(knd) x, cqsqrt
        cqsqrt = sqrt(x)
        return
    end function cqsqrt

    function iqint(x)
        real(knd) x
        integer iqint
        iqint = int(x)
        return
    end function iqint

    function qreal(x)
        complex(knd) x
        real(knd) qreal
        qreal = real(x)
        return
    end function qreal

    function qimag(x)
        complex(knd) x
        real(knd) qimag
        qimag = imag(x)
        return
    end function qimag

    function qcmplx(x, y)
        complex(knd) qcmplx
        real(knd) x, y
        qcmplx = cmplx(x, y, knd)
        return
    end function qcmplx

    function qint(x)
        real(knd) x, qint
        qint = aint(x, knd)
    end function qint

end module elemfunc