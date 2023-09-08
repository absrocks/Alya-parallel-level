/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#ifndef DISTRIBUTED_VOLUME_APPLICATION_H
#define DISTRIBUTED_VOLUME_APPLICATION_H

#define NOISE_SIZE  64

//----------------------------------------------------------------------
/// Distributed volume application class
/// A sample volume manager class that creates an analytic volume (it doesn't allocate data)
/// a split it in bricks of equal size to be used in data distributed environments.
class Distributed_volume_application
{    
public:

    enum Volume_type
    {
        DEFAULT     = 0,
        TURBULENCE,
    };

//----------------------------------------------------------------------
/// Construct the distributed volume
/// \param[in] volume_size  The target volume size
/// \param[in] brick_size   The brick size used for distribution
///
    Distributed_volume_application(
        const mi::math::Vector<mi::Uint32, 3>   &volume_size, 
        const mi::math::Vector<mi::Uint32, 3>   &brick_size,
        Volume_type                             volume_type = DEFAULT);
    
//----------------------------------------------------------------------
/// Get the total number of bricks in the volume
/// \return     The number of bricks
///
    mi::Uint32 nb_bricks() const;
    
//----------------------------------------------------------------------
/// Get a volume brick given its 1D brick index
/// \param[in] brick_idx    The brick index (0 to nb_bricks() - 1)
/// \param[out] bbox        The brick bounding box
/// \param[out] brick       The destination brick buffer
///
    void get_brick(
        mi::Uint32                      brick_idx, 
        mi::math::Bbox<mi::Sint32, 3>   &bbox, 
        mi::Uint8                       *brick) const;
        
    void get_brick(
        mi::Uint32                      brick_idx, 
        mi::math::Bbox<mi::Float32, 3>   &bbox, 
        mi::Uint8                       *brick) const;
        
        
private:

    void generate_noise();
    
    mi::Float32 smooth_noise(
        mi::Float32 x, 
        mi::Float32 y, 
        mi::Float32 z) const;
        
    mi::Float32 turbulence(
        mi::Float32 x, 
        mi::Float32 y, 
        mi::Float32 z, 
        mi::Float32 size) const;

    mi::math::Vector<mi::Uint32, 3>   m_volume_size;
    mi::math::Vector<mi::Uint32, 3>   m_brick_size;
    mi::math::Vector<mi::Uint32, 3>   m_volume_bricks_size;
    Volume_type                       m_volume_type;
    mi::Float32                       m_noise[NOISE_SIZE][NOISE_SIZE][NOISE_SIZE];
};

#endif // DISTRIBUTED_VOLUME_APPLICATION_H

