program simanhosc
    use, intrinsic :: ieee_arithmetic
    implicit none
    
    character(256) :: Foldername
    real(4) :: m, new_m, varphi, partition
    integer :: k

    m = 50.0
    partition = m / 25

    Foldername = 'Sch5' 
    call Create_Folder(Foldername)


    Foldername = 'Sch5/BH0'
    Varphi = 0
    call Simulate_BH(m, varphi, Foldername, partition)

    do k = 1, 10, 1
        varphi = k  
        write(Foldername, '(A,I0)') 'Sch5/BH', k
        new_m = Find_Q_from_M_fixed_Hayward(m, varphi)

        call Simulate_BH(new_m, varphi, Foldername, partition)

    end do

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! implicit Gauss-Legendre methods; symplectic with arbitrary Hamiltonian, A-stable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! equations of motion General Kerr space-time
    subroutine evalf(y, dydx, t, q, r0)
        real y(3), dydx(3), t, q,r0                     ! Do I change it?
        ! y(1) = t(r), y(2) = phi(r)
        if (t*t*t - q*q*t + q*q*r0 > 0) then
            dydx(1) = +t*t*t/((t-r0)*sqrt((t*(t*t*t - q*q*t + q*q*r0))) ) !t(r) r intercambiado por t
            dydx(2) = -q/sqrt((t*(t*t*t - q*q*t + q*q*r0) )) 
            dydx(3) = 0
        elseif (t*t*t - q*q*t + q*q*r0 < 0) then
            dydx(1) = -t*t*t/((t-r0)*sqrt(abs(t*(t*t*t - q*q*t + q*q*r0))) ) !t(r) r intercambiado por t
            dydx(2) = +q/sqrt(abs(t*(t*t*t - q*q*t + q*q*r0) )) 
            dydx(3) = 0
        end if
    end subroutine evalf

    ! 8th order implicit Gauss-Legendre integrator
    ! only change n -> the number of equations
    subroutine gl8(y, dt,t, q, r0)
        integer, parameter :: s = 4, n = 3
        real y(n), g(n,s), dt,  t, q, r0; integer i, k, j
        
        ! Butcher tableau for 8th order Gauss-Legendre method
        real, parameter :: a(s,s) = reshape((/ &
                 0.869637112843634643432659873054998518Q-1, -0.266041800849987933133851304769531093Q-1, &
                 0.126274626894047245150568805746180936Q-1, -0.355514968579568315691098184956958860Q-2, &
                 0.188118117499868071650685545087171160Q0,   0.163036288715636535656734012694500148Q0,  &
                -0.278804286024708952241511064189974107Q-1,  0.673550059453815551539866908570375889Q-2, &
                 0.167191921974188773171133305525295945Q0,   0.353953006033743966537619131807997707Q0,  &
                 0.163036288715636535656734012694500148Q0,  -0.141906949311411429641535704761714564Q-1, &
                 0.177482572254522611843442956460569292Q0,   0.313445114741868346798411144814382203Q0,  &
                 0.352676757516271864626853155865953406Q0,   0.869637112843634643432659873054998518Q-1 /), (/s,s/))
        real, parameter ::   b(s) = (/ &
                 0.173927422568726928686531974610999704Q0,   0.326072577431273071313468025389000296Q0,  &
                 0.326072577431273071313468025389000296Q0,   0.173927422568726928686531974610999704Q0  /)
        
        ! iterate trial steps
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dt, g(:,i), t, q, r0)
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
        if (t*t*t - q*q*t + q*q*r0 > 0) then
            t = t + dt
        end if
        if (t*t*t - q*q*t + q*q*r0 < 0) then
            t = t -dt
        end if
    end subroutine gl8

    function Final_cond(t_initial, q, r0) 
        real :: Final_cond(3), t_initial, q, r0

        Final_cond(1) = t_initial + r0 * log(t_initial/r0)

        Final_cond(2) = q / t_initial

        Final_cond(3) = 0
    end function Final_cond

    !!!! Functions and routines regarding the manage of the data/
    subroutine Create_Folder(Folder)
        implicit none
        character(256), intent(in) :: Folder
        logical :: folder_exists
        integer :: status
        character(256) :: command
    
        ! Verificar si la carpeta ya existe
        inquire(file=trim(Folder), exist=folder_exists)
    
        if (folder_exists) then
            ! Eliminar la carpeta y su contenido en Linux
            command = 'rm -rf "' // trim(Folder) // '"'
            status = system(command)
            if (status /= 0) then
                print *, "Error al eliminar la carpeta: ", trim(Folder)
            endif
        endif
    
        ! Crear la carpeta en Linux
        command = 'mkdir "' // trim(Folder) // '"'
        status = system(command)
        if (status /= 0) then
            print *, "Error al crear la carpeta: ", trim(Folder)
        endif
    end subroutine Create_Folder

    subroutine read_file(File, Folder, max_rows, n, df, no_cols)
        implicit none
        character(256), intent(in) :: file, folder
        character(256) :: file_temp
        integer, intent(in) :: max_rows, no_cols
        integer, intent(out) :: n
        real(4), allocatable, intent(out) :: df(:, :)
        
        integer :: status, iostat_code, unit_id, j

        unit_id = 10

        ! Asignar espacio para los arrays dinámicos de los datos originales
        allocate(df(max_rows, no_cols))

        write(file_temp, '(A,"/",A,".txt")') trim(folder), trim(file)
        open(unit=unit_id, file=trim(file_temp), status='unknown', iostat=iostat_code, action="read")
        if (iostat_code /= 0) then
            print *, "Error al abrir el archivo:", trim(file_temp)
            print *, "Código de error IOSTAT:", iostat_code
            return
        endif

        ! Leer los datos en los arrays
        n = 0
        do
            read(unit_id, *, iostat=status) (df(n+1, j), j=1, no_cols)
            if (status /= 0) exit
            n = n + 1
        end do
        close(unit_id)
    end subroutine read_file
    
    subroutine fill_with_intermediate_points(File, Folder, df, n, dt_min, threshold, no_cols)
        implicit none
        character(256), intent(in) :: File, Folder
        real(4), intent(in) :: threshold, dt_min
        integer, intent(inout) :: n
        real, intent(inout) :: df(:,:)
        integer, intent(in) :: no_cols

        character(256) :: file_temp
        integer :: i, j, max_rows, k, l, m, iostat_code, unit_id
        real(4), allocatable :: new_df(:,:), diff(:), var_dt(:)

        unit_id = 20

        max_rows = 100000
        allocate(new_df(max_rows, no_cols), diff(no_cols), var_dt(no_cols))

        write(file_temp, '(A,"/",A , ".txt")') trim(folder), trim(file)
        
        open(unit=unit_id, file=trim(file_temp), status="replace", action="write", iostat=iostat_code)
        if (iostat_code /= 0) then
            print *, "Error al abrir el archivo:", trim(file_temp)
            print *, "Código de error IOSTAT:", iostat_code
            return
        endif

        j = 0  ! Índice para la nueva lista con valores intermedios
        ! Procesar los datos para agregar valores intermedios
        do i = 1 , n - 1
            ! Copiar el valor original
            j = j + 1
            do l = 1, no_cols, 1
                new_df(j, l) = df(i, l)
            end do

            ! Evaluar si la diferencia en la primera columna excede el umbral
            if (abs(df(i+1, 4) - df(i, 4)) > threshold) then
                m =  ceiling(abs(df(i+1, 4) - df(i, 4))/dt_min)
                do l = 1, no_cols, 1
                    diff(l) = df(i+1, l) - df(i, l)
                    var_dt(l) = diff(l) / m
                end do

                do k = 1, m-1
                    j = j + 1
                    do l = 1, no_cols, 1
                        new_df(j, l) = df(i, l) + k * var_dt(l)
                    end do
                end do    
            endif
        end do

        ! Copiar el último punto
        j = j + 1
        do l = 1, no_cols, 1
            new_df(j, l) = df(n, l)
        end do
        n = j

        do i = 1, n
            write(unit_id, *) (new_df(i, j), j = 1, no_cols)
        end do

        df = new_df
        ! Liberar memoria de los arrays 
        deallocate(new_df, diff, var_dt)
        close(unit_id)
    end subroutine fill_with_intermediate_points
    
    subroutine shrink_data(File, Folder, which, df, n, min, max, no_cols)
        implicit none
        character(256), intent(in) :: File, Folder
        character(1), intent(in):: which
        real(4), intent(in) :: min, max
        real(4), intent(inout) :: df(:, :)
        integer, intent(inout) :: n
        integer, intent(in) :: no_cols

        character(256) :: file_temp
        integer :: unit_id, i, j, n_temp, column
        real(4), allocatable :: df_temp(:, :)

        unit_id = 30

        write(file_temp, '(A,"/",A , ".txt")') trim(folder), trim(file)
        ! print *, "Nombre del archivo de Salida: ", trim(file_temp)

        ! call system('del ' // trim(file_temp))  !tehre are some problems with the use of system
        call execute_command_line('rm "' // trim(file_temp) // '"')
        ! print*, df(1,1)
        ! print*, "n", n, no_cols
        allocate(df_temp(n, no_cols))

        select case (which)
            case ('X')
                column = 1
            case ('Y')
                column = 2
            case ('Z')
                column = 3
            case ('T')
                column = 4
            case default
                print *, "Error: 'which' must be one of 'X', 'Y', 'Z', or 'T'."
                stop
        end select

        n_temp = 0
        do i = 1, n
            if (df(i, column) >= min .and. df(i, column) <= max) then
                n_temp = n_temp + 1
                do j = 1, no_cols
                    df_temp(n_temp, j) = df(i, j)
                end do
            endif
        end do

        n = n_temp

        do i = 1, n_temp, 1
            do j = 1, no_cols, 1
                df(i, j) = df_temp(i, j)
            end do
        end do
        deallocate(df_temp)

        open(unit=unit_id, file=trim(file_temp), status="replace")

        if (n == 0) then
            call execute_command_line('rm "' // trim(file_temp) // '"')
        else 
            do i = 1, n
                write(unit_id, *) (df(i, j), j = 1, no_cols)
            end do
        end if

        close(unit_id)
    end subroutine shrink_data

    subroutine Backwards_Integration(Folder, varphi_start, varphi_end, varphi_step, P_start, P_end, P_step, m, t_initial, t_final, dt)
        character(256), intent(in) :: Folder
        real(4), intent(in) :: varphi_end, varphi_start, varphi_step
        real(4), intent(in) :: P_end, P_start, P_step
        real(4), intent(in) :: m, t_initial, t_final, dt

        character(256) :: foldername_varphi, full_subfolder, filename
        integer :: varphi_value, P_value, j
        real(4) :: P_alpha, varphi, dt1, var_m, x, y1, z, t, y(3), y_previous(3), tol_x, t_previous, P_step_var
        
        do varphi_value = 0, floor((varphi_end - varphi_start)/varphi_step) - 1,  1
            varphi = varphi_start + varphi_value * varphi_step

            write(foldername_varphi, '(A,I0)') 'phi', int(varphi_value)
            full_subfolder = trim(Folder) // '/' // trim(foldername_varphi)
            call Create_Folder(full_subfolder)

            ! do P_value = 0, floor((P_end - P_start)/P_step), 1
                ! P_alpha = P_start + P_value * P_step
            P_value = 0
            P_alpha = P_start
            do
                P_value = P_value + 1
                if ((P_alpha < P_start) .or. (P_alpha > P_end)) exit
                if (abs(abs(P_alpha) - 2.22864* (2*m)) < 0.4 * (2*m)) then
                    P_alpha = P_alpha + P_step * (0.1 + (abs(abs(P_alpha) - 2.22864* (2*m)) * (0.9/(0.4 * 2*m))))
                else
                    if (P_value == 1) then
                        P_alpha = P_alpha
                    else
                        P_alpha = P_alpha + P_step
                    end if
                end if

                ! if (P_value == floor((P_end - P_start)/P_step) + 1) then
                !     P_alpha = -2.22864 * 2*m
                ! else if (P_value == floor((P_end - P_start)/P_step) + 2) then
                !     P_alpha = 2.22864 * 2*m
                ! endif

                ! Generar el nombre del archivo dinámicamente 
                write(filename, '(A, "/data", F8.2, ".txt")') trim(full_subfolder), P_alpha
                ! Abrir el archivo para escribir los datos
                open(10, file=trim(filename), status="replace")
                
                y = Final_cond(t_initial, P_alpha , 2*m)
                var_m = m
                t = t_initial

                ! Bucle para calcular los valores de x, z, y t
                do j = 1, 10000
                    y_previous = y
                    t_previous = t

                    if (j == 1) then
                        x = t * sin(y(2))
                        z = t * cos(y(2))
                        write(10,*) x, 0, z, y(1), y(2), t, dt1, var_m
                    end if

                    if (t >= m * 2.5) then
                        dt1 = -dt / 5
                    else if (t < m * 2.5) then
                        if (y(2) > 1.5707) then
                            ! dt1= -dt / (50 + 700 * (y(2) + 1.5708) / 1.5708)   
                            dt1 = -dt/ ((10*m - t)) 
                        else if (y(2) <= 1.5707) then
                            dt1 = -dt/ ((10*m - t)) 
                        end if      
                    end if

                    call gl8(y, dt1, t, P_alpha, 2*var_m)   
                    x = t * sin(y(2))
                    z = t * cos(y(2))
                    
                    if (P_alpha < 0 .and. x > 0) then
                        y = y_previous
                        x = t * sin(y(2))
                        z = t * cos(y(2))
                        t = t_previous
                        do  
                            dt1 = -dt / (200 + (2.5 - x)*(200/2.5))
                            call gl8(y, dt1, t, P_alpha, 2*var_m)   
                            x = t * sin(y(2))
                            z = t * cos(y(2))
                            if (x > 0) exit
                            write(10,*) x, 0, z, y(1), y(2), t, dt1, var_m
                        end do
                    else if (P_alpha > 0 .and. x < 0) then
                        y = y_previous
                        x = t * sin(y(2))
                        z = t * cos(y(2))
                        t = t_previous
                        do  
                            dt1 = -dt / 250
                            call gl8(y, dt1, t, P_alpha, 2*var_m)   
                            x = t * sin(y(2))
                            z = t * cos(y(2))
                            if (x < 0) exit
                            write(10,*) x, 0, z, y(1), y(2), t, dt1, var_m
                        end do
                    end if

                    if ((P_alpha > 0 .and. x < 0) .or. &
                        (P_alpha < 0 .and. x > 0) .or. &
                        (y(1) < t_final) .or. &
                        ieee_is_nan(y(1)) .or. &
                        ieee_is_nan(y(2)) .or. &
                        ieee_is_nan(y(3)) .or. &
                        (y(1) > y_previous(1))) then
                        exit
                    end if

                    write(10,*) x, 0, z, y(1), y(2), t, dt1, var_m
                    var_m = m * t**3 / (t**3 + (varphi)**3)

                end do
                
                close(10)  ! Cerrar el archivo después de cada valor de q
                print *, 'Calculos completados para P_alpha =', P_alpha, "q = ", varphi, "Masa = ", m
            end do 
        end do
    end subroutine Backwards_Integration

    subroutine Improve_data(Folder, varphi_start, varphi_end, varphi_step, P_start, P_end, P_step, m, t_initial, t_final, dt, dt_min, threshold)
        character(256), intent(in) :: Folder
        real(4), intent(in) :: varphi_end, varphi_start, varphi_step
        real(4), intent(in) :: P_end, P_start, P_step
        real(4), intent(in) :: m, t_initial, t_final, dt, dt_min, threshold

        character(256) :: foldername_varphi, full_subfolder, filename
        real(4), allocatable :: df(:, :)
        real(4) :: P_alpha, varphi, min, max, P_step_var
        integer :: no_cols, max_rows, n, varphi_value, P_value

        no_cols = 8

        do varphi_value = 0, floor(abs(varphi_end - varphi_start)/varphi_step) - 1, 1
            varphi = varphi_start + varphi_value * varphi_step
            ! do P_value = 0, floor((P_end - P_start)/P_step) , 1
                ! P_alpha = P_start + P_value * P_step
            P_value = 0
            P_alpha = P_start
            do
                P_value = P_value + 1
                if ((P_alpha < P_start) .or. (P_alpha > P_end)) exit
                if (abs(abs(P_alpha) - 2.22864* (2*m)) < 0.4 * (2*m)) then
                    P_alpha = P_alpha + P_step * (0.1 + (abs(abs(P_alpha) - 2.22864* (2*m)) * (0.9/(0.4 * 2*m))))
                else
                    if (P_value == 1) then
                        P_alpha = P_alpha
                    else
                        P_alpha = P_alpha + P_step
                    end if
                end if
                
                ! if (P_value == floor((P_end - P_start)/P_step) + 1) then
                !     P_alpha = -2.22864 * 2*m
                ! else if (P_value == floor((P_end - P_start)/P_step) + 2) then
                !     P_alpha = 2.22864 * 2*m
                ! endif

                write(foldername_varphi, '(A,I0)') 'phi', int(varphi_value)
                full_subfolder = trim(Folder) // '/' // trim(foldername_varphi)
                max_rows = 100000

                write(filename, '("data", F8.2)') P_alpha

                call read_file(filename, full_subfolder, max_rows, n, df, no_cols)
                call fill_with_intermediate_points(filename, full_subfolder, df, n, dt_min, threshold, no_cols)

                min = -m*20 - 10
                max = m*30+ 10
                call shrink_data(filename, full_subfolder, 'T', df, n, min, max, no_cols)
                min = -m*25 - 10
                max = m*25 + 10
                call shrink_data(filename, full_subfolder, 'X', df, n, min, max, no_cols)
                print*, 'Datos mejorados para P_alpha =', P_alpha, "Varphi = ", varphi, "Masa = ", m
            end do
    end do
    end subroutine Improve_data

    subroutine Simulate_BH(m, varphi, Foldername, partition)
        real(4), intent(in) :: m, varphi, partition
        character(256), intent(in) ::  Foldername

        ! For integration
        integer P_value, varphi_value
        real(4) :: Lz, E, alpha, dt, t_final, t_initial
        real(4) :: P_end, P_start, P_step
        real(4) :: varphi_start, varphi_end, varphi_step
        real(4) :: max_iter, tol_x

        ! For improving the data
        real(4) :: dt_min, threshold, minX, maxX, minZ, maxZ, minT, maxT, q_bbh_sbh

        ! For interaction
        integer :: user_choice

        !! Change Here 
        E = 1
        t_initial = m * 100
        t_final = - m * 30
        dt = m / 10
        dt_min = m / 100
        threshold = dt_min

        q_bbh_sbh = 2.67848 * (2 * m)
        ! q_bbh_sbh = 10*m

        P_end = q_bbh_sbh
        P_start = -P_end
        P_step = partition
        if (varphi == 0) then
            varphi_start = 0
            varphi_end = 1
            varphi_step = 1
        else 
            varphi_start = varphi
            varphi_end = 2 * varphi
            varphi_step =  varphi
        end if

        call Create_Folder(foldername)
        call Backwards_Integration(foldername, varphi_start, varphi_end, varphi_step, P_start, P_end, P_step, m, t_initial, t_final, dt)        
        call Improve_data(foldername, varphi_start, varphi_end, varphi_step, P_start, P_end, P_step, m, t_initial, t_final, dt, dt_min, threshold)   
    end subroutine Simulate_BH

    function Find_Q_from_M_fixed_Hayward(m, Varphi)
        real(4) :: Find_Q_from_M_fixed_Hayward, m, varphi
        
        Find_Q_from_M_fixed_Hayward = m * ((2 * m)**3 + varphi**3) / (2 * m)**3
    end function Find_Q_from_M_fixed_Hayward


end

    
