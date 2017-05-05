        program variableTesting
            implicit none
            ! declare variables
            integer :: total
            real :: average
            complex :: cx
            logical :: done
            character(len=80) :: message !string of 80 characters
        
            ! assigning values
            total = 20000
            average = 1666.67
            done = .true.
            message = "Hello World!"
            cx = (3.0, 5.0) ! complex 3.0 + 5.0i
            
            print *, total
            print *, average
            print *, cx
            print *, done
            print *, message
            
        end program variableTesting