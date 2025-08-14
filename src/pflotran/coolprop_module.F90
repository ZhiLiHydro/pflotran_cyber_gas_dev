MODULE CoolProp_module
    INTERFACE
        FUNCTION PropsSI (output, name1, prop1, name2, prop2, fluidname) BIND(C, NAME='PropsSI')
            use iso_c_binding
            real(C_DOUBLE) :: PropsSI
            character(KIND=c_char) :: output(*)
            character(c_char) :: name1(*)
            real(C_DOUBLE), VALUE :: prop1
            character(c_char) :: name2(*)
            real(C_DOUBLE), VALUE :: prop2
            character(kind=c_char) :: fluidname(*)
        END FUNCTION PropsSI
        FUNCTION HAPropsSI (output, name1, prop1, name2, prop2, name3, prop3) BIND(C, NAME='HAPropsSI')
            use iso_c_binding
            real(C_DOUBLE) :: HAPropsSI
            character(KIND=c_char) :: output(*)
            character(c_char) :: name1(*)
            real(C_DOUBLE), VALUE :: prop1
            character(c_char) :: name2(*)
            real(C_DOUBLE), VALUE :: prop2
            character(c_char) :: name3(*)
            real(C_DOUBLE), VALUE :: prop3
        END FUNCTION HAPropsSI
                FUNCTION get_input_pair_index(input) BIND(C, NAME='get_input_pair_index')
            use iso_c_binding
            integer(C_long) :: get_input_pair_index
            character(KIND=c_char) :: input(*)
            
        END FUNCTION get_input_pair_index
        
        FUNCTION get_param_index(input) BIND(C, NAME='get_param_index')
            use iso_c_binding
            integer(C_long) :: get_param_index
            character(KIND=c_char) :: input(*)
            
        END FUNCTION get_param_index
        
        FUNCTION AbstractState_factory(backend, fluidname, errcode, buffer, buffersize) BIND(C, NAME='AbstractState_factory')
            use iso_c_binding
            integer(C_long) :: AbstractState_factory
            character(kind=c_char):: backend(*)
            character(kind=c_char) :: fluidname(*)
            integer(C_long) :: errcode
            character(kind=c_char) :: buffer(*)
            integer(C_long), VALUE :: buffersize
            
        END FUNCTION AbstractState_factory
        
        SUBROUTINE AbstractState_update(handle, input, prop1, prop2, errcode, &
            buffer, buffersize) BIND(C, NAME='AbstractState_update')
            use iso_c_binding
            integer(C_long), VALUE :: handle
            integer(C_long), VALUE :: input
            real(C_DOUBLE), VALUE :: prop1
            real(C_DOUBLE), VALUE :: prop2
            integer(C_long) :: errcode
            character(kind=c_char) :: buffer(*)
            integer(C_long), VALUE :: buffersize
            
        END SUBROUTINE AbstractState_update
        
        FUNCTION AbstractState_keyed_output(handle, output, errcode, &
            buffer, buffersize) BIND(C, NAME='AbstractState_keyed_output')
            use iso_c_binding
            real(C_DOUBLE) :: AbstractState_keyed_output
            integer(C_long),VALUE :: handle
            integer(C_long),VALUE :: output
            integer(C_long) :: errcode
            character(kind=c_char) :: buffer(*)
            integer(C_long), VALUE :: buffersize
            
        END FUNCTION AbstractState_keyed_output
        
        FUNCTION AbstractState_first_partial_deriv(handle, Of, Wrt, Constnt, errcode, &
            buffer, buffersize) BIND(C, NAME='AbstractState_first_partial_deriv')
            use iso_c_binding
            real(C_DOUBLE) :: AbstractState_first_partial_deriv
            integer(C_long),VALUE :: handle
            integer(C_long),VALUE :: Of
            integer(C_long),VALUE :: Wrt
            integer(C_long),VALUE :: Constnt
            integer(C_long) :: errcode
            character(kind=c_char) :: buffer(*)
            integer(C_long), VALUE :: buffersize

        END FUNCTION AbstractState_first_partial_deriv

    END INTERFACE
END MODULE CoolProp_module