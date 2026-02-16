! This module defines variable precision for all common data types.
module kinds
    implicit none

    ! variable types
    integer, parameter :: byte = selected_int_kind(0)
    integer, parameter :: i4   = selected_int_kind(6) 
    integer, parameter :: long = selected_int_kind(18)
    integer, parameter :: r4   = selected_real_kind(6)
    integer, parameter :: r8   = selected_real_kind(13)

    ! zero-values
    real(r8), parameter      :: dzero = 0.0
    real(r4), parameter      :: rzero = 0.0
    integer(i4), parameter   :: izero = 0
    integer(byte), parameter :: bzero = 0

    ! one-values
    real(r8), parameter      :: done = 1.0
    real(r4), parameter      :: rone = 1.0
    integer(i4), parameter   :: ione = 1
    integer(byte), parameter :: bone = 1

    ! misc-values
    integer(byte), parameter :: btwo = 2
    integer(byte), parameter :: bnegone  = -1
    integer(byte), parameter :: bnegfive = -5

end module kinds
