PROGRAM main
    implicit none
    double precision, dimension (:, :), allocatable :: A
    double precision, dimension (:), allocatable :: tau, work, exmp_work
    integer :: n, i, j, info, lda, ihi, ilo, lwork
    real ::  start, end

    n = 100

    allocate (A(n, n), tau(n), exmp_work(n))

    do j = 1, n
        do i = 1, n
            A(i, j) = i + j
        end do
    end do
    lda = n
    lwork = -1
    ilo = 1
    ihi = n

    call dgehrd(n, 1, n, A, n, tau, exmp_work, lwork, info)
    lwork = exmp_work(1)

    allocate(work(lwork))
    call cpu_time(start)

    call dgehrd(n, 1, n, A, n, tau, work, lwork, info)

    call cpu_time(end)
    write (*,*) "Time in seconds: "
    print *, end - start

    deallocate (A, tau, work, exmp_work)

END PROGRAM main
