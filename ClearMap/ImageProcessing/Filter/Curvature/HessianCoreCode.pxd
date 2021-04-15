from numpy cimport uint8_t, uint16_t, float32_t, double_t


ctypedef fused dtype_t_source:
    uint8_t
    uint16_t
    float32_t
    double_t


ctypedef fused dtype_t_sink:
    uint8_t
    uint16_t
    float32_t
    double_t


cdef void core(void kernel(dtype_t_sink*, 
                           double, double, double, double) nogil,
               dtype_t_source[:, :, ::1]  source,
               dtype_t_sink[:, :, :, ::1] sink,
               double par) except *