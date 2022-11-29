program main
    implicit none
    double precision, dimension(:,:), allocatable :: A
    integer :: i, j

    integer, parameter :: n = 4

    allocate (A(n, n))
    !doubleprecision, dimension(n, n) :: A
    !data A /4, 2, 3, 1, 2, -3, 5, 1, 3, 1, 2, 1, 3, 1, 1, 2/
    call random_number(A)

    call hausholders_transformation(A, n)
    !write(*, *) A
    do i=1,n
        write(*,*)((anint(A(j,i)*100)/100),j=1,n)
    end do
    deallocate (A)
END PROGRAM main

subroutine hausholders_transformation(A, n)
    implicit none
    integer :: i, j, t
    integer, parameter :: k = 2
    integer, intent(in) :: n
    DOUBLE PRECISION, intent(inout), dimension(n, n) :: A
    double precision s, norm, norm_x, b
    double precision, dimension(:), allocatable :: x, y, z
    double precision, dimension(:, :), allocatable :: tmp, tmp1, prom
    allocate (tmp(n, n), tmp1(n, n), prom(n, n), x(n), y(n), z(n))

    do i = 1, n - 2
        s = 0
        norm = 0
        norm_x = 0
        b = 0

        do j = i + 2, n
            s = s + (A(i, j) ** 2)
        end do

        norm = sqrt(abs(A(i, i+1))** 2 + s)
        do j = 1, n
            if (j < i + 1) then
                x(j) = 0
            else if (j == i + 1) then
                x(j) = A(i, j) - norm
            else
                x(j) = A(i, j)
            end if
        end do
        
        norm_x = sqrt(abs(x(i + 1))** 2 + s)
        do j = 1, n
            x(j) = x(j) / norm_x
        end do

        call dot_matr_vec1(A, x, y, n)
        call dot_matr_vec2(A, x, z, n)

        call dot(x, y, tmp, n)
        call dot(z, x, tmp1, n)

        do j = 0, n
            b = b + y(j) * x(j)
        end do
        b = b * 4

        call dot1(x, x, prom, n ,b)

        do j = 1, n
            do t = 1, n
                A(t, j) = A(t, j) - tmp1(t, j) - tmp(t, j) + prom(t, j)
            end do
        end do

    end do

end subroutine hausholders_transformation

subroutine dot_matr_vec1(a, b, ans, n)
    implicit none
    integer :: i, j
    integer, intent(in) :: n
    DOUBLE PRECISION, INTENT(IN), dimension(n) :: b
    DOUBLE PRECISION, INTENT(IN), dimension(n, n) :: a

    DOUBLE PRECISION, INTENT(INOUT), dimension(n) :: ans
    do i = 1, n
        ans(i) = 0
        do j = 1, n
            ans(i) = ans(i) + a(i, j) * b(j)
        end do
    end do
end subroutine dot_matr_vec1


subroutine dot_matr_vec2(a, b, ans, n)
    implicit none
    integer :: i, j
    integer, intent(in) :: n
    DOUBLE PRECISION, INTENT(IN), dimension(n, n) :: a
    DOUBLE PRECISION, INTENT(IN), dimension(n) :: b

    DOUBLE PRECISION, INTENT(INOUT), dimension(n) :: ans
    do i = 1, n
        ans(i) = 0
        do j = 1, n
            ans(i) = ans(i) + a(j, i) * b(j)
        end do
    end do
end subroutine dot_matr_vec2


subroutine dot(a, b, ans, n)
    implicit none
    integer :: i, j
    integer, intent(in) :: n
    DOUBLE PRECISION, INTENT(IN), dimension(n) :: a, b
    DOUBLE PRECISION, INTENT(INOUT), dimension(n, n) :: ans
    do i = 1, n
        do j = 1, n
            ans(j, i) = 2 * a(i) * b(j)
        end do
    end do
end subroutine dot

subroutine dot1(a, b, ans, n, k)
    implicit none
    integer :: i, j
    integer, intent(in) :: n
    double precision, intent(in) :: k
    DOUBLE PRECISION, INTENT(IN), dimension(n) :: a, b
    DOUBLE PRECISION, INTENT(INOUT), dimension(n, n) :: ans
    do i = 1, n
        do j = 1, n
            ans(j, i) = k * a(i) * b(j)
        end do
    end do
end subroutine dot1