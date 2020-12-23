! Created by drakosha on 23.12.2020.

program test_inverse
    use matrix
    use regime
    implicit none

    complex(knd), allocatable, dimension(:,:) :: A, B
    integer :: n, i, j

    open(10, file = "../testing/tests/inverse_matrix_test.txt")
    read(10, *) n
    write(*,*) 'n = ', n
    allocate(A(n,n), B(n,n))
    read(10, *) A
    call inverse_matrix(A, n, B)
    write(*,*) B
    call quick_inverse_matrix(A, n, B)
    write(*,*) B

end program test_inverse