    subroutine MCMCPardiso(nnode, nelmt,  enode,nbc, id, outout)
    implicit none
    integer nnode, nelmt

    ! element data
    integer enode(nelmt,0:9), etype(nelmt)

    ! boundary condition data
    integer nbc(3), id(nnode,2), ndof, outout

    ! frequency data


    integer itmp, jtmp, inode, jnode, idof, jdof, ielmt
    integer,allocatable :: ia(:), ja(:), idcheck(:)
    integer equation, node, position



    allocate (ia(nbc(2)+1),idcheck(nbc(2)))

    ia(1) = 1

    do equation = 1 , nbc(2)

        node = 0

        do inode = 1 , nnode
            do idof = 1 , 2

                if (id(inode,idof) == equation+nbc(1)) then

                    node = inode

                endif

            enddo
        enddo

        idcheck(:) = 0

        do ielmt = 1 , nelmt
            do inode = 1 , enode(ielmt,0)

                if (enode(ielmt,inode) == node) then

                    do jnode = 1 , enode(ielmt,0)

                        do jdof = 1 , 2

                            if ((id(enode(ielmt,jnode),jdof) >= nbc(1)+1) .and. &
                                (id(enode(ielmt,jnode),jdof) <= nbc(1)+nbc(2))) then

                                idcheck(id(enode(ielmt,jnode),jdof)-nbc(1)) = 1

                            endif

                        enddo

                    enddo

                endif

            enddo
        enddo

        ia(equation+1) = ia(equation)

        do idof = 1 , nbc(2)

            if (idcheck(idof) == 1) then

                ia(equation+1) = ia(equation+1) + 1

            endif

        enddo

    enddo

    allocate (ja(ia(nbc(2)+1)-1))

    position = 0

    do equation = 1 , nbc(2)

        node = 0

        do inode = 1 , nnode
            do idof = 1 , 2

                if (id(inode,idof) == equation+nbc(1)) then

                    node = inode

                endif

            enddo
        enddo

        idcheck(:) = 0

        do ielmt = 1 , nelmt
            do inode = 1 , enode(ielmt,0)

                if (enode(ielmt,inode) == node) then

                    do jnode = 1 , enode(ielmt,0)

                        do jdof = 1 , 2

                            if ((id(enode(ielmt,jnode),jdof) >= nbc(1)+1) .and. &
                                (id(enode(ielmt,jnode),jdof) <= nbc(1)+nbc(2))) then

                                idcheck(id(enode(ielmt,jnode),jdof)-nbc(1)) = 1

                            endif

                        enddo

                    enddo

                endif

            enddo
        enddo

        do idof = 1 , nbc(2)

            if (idcheck(idof) == 1) then

                position = position + 1

                ja(position) = idof

            endif

        enddo

    enddo
    do itmp = 1 , nbc(2)+1

        write (outout,*) ia(itmp)

    enddo

    write (outout,*)

    do itmp = 1 , ia(nbc(2)+1)-1

        write (outout,*) ja(itmp)

    enddo

    deallocate (ia, ja)

    end subroutine MCMCPardiso