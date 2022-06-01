/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include <mi/dice.h>

#include "distributed_volume_application.h"

#include "common/forwarding_logger.h"

//----------------------------------------------------------------------
Distributed_volume_application::Distributed_volume_application(
    const mi::math::Vector<mi::Uint32, 3>   &volume_size, 
    const mi::math::Vector<mi::Uint32, 3>   &brick_size,
    Volume_type                             volume_type)
{
    m_volume_size = volume_size;
    m_brick_size = brick_size;

    m_volume_bricks_size.x = volume_size.x / brick_size.x + 
        ((volume_size.x%brick_size.x != 0) ? 1 : 0);

    m_volume_bricks_size.y = volume_size.y / brick_size.y + 
        ((volume_size.y%brick_size.y != 0) ? 1 : 0);

    m_volume_bricks_size.z = volume_size.z / brick_size.z + 
        ((volume_size.z%brick_size.z != 0) ? 1 : 0);
        
    m_volume_type = volume_type;
    
    if(volume_type == TURBULENCE)
        generate_noise();
}
        
//----------------------------------------------------------------------
mi::Uint32 Distributed_volume_application::nb_bricks() const
{
    return m_volume_bricks_size.x*m_volume_bricks_size.y*m_volume_bricks_size.z;
}

//----------------------------------------------------------------------
void Distributed_volume_application::get_brick(
    mi::Uint32                      brick_idx, 
    mi::math::Bbox<mi::Sint32, 3>   &bbox, 
    mi::Uint8                       *brick) const
{
    mi::Uint32 idx = brick_idx;
    // convert linear brick index to x,y,z position
    // get z
    mi::Uint32 z = idx / (m_volume_bricks_size.x*m_volume_bricks_size.y);
    idx = idx - z*(m_volume_bricks_size.x*m_volume_bricks_size.y);
    
    // get y
    mi::Uint32 y = idx / m_volume_bricks_size.x;

    // get x
    mi::Uint32 x = idx - y*m_volume_bricks_size.x;
    
    //calculate brick bounding box
    bbox.min.x = x*m_brick_size.x;
    bbox.min.y = y*m_brick_size.y;
    bbox.min.z = z*m_brick_size.z;
    
    bbox.max.x = bbox.min.x + m_brick_size.x;
    bbox.max.y = bbox.min.y + m_brick_size.y;
    bbox.max.z = bbox.min.z + m_brick_size.z;
    
    // INFO_LOG << "Brick (" << brick_idx << ") coordinates: " << x << ", " << y << ", " << z <<
        // " and BBox: " << bbox;
        
    // calculate data
    if(brick != NULL)
    {
        switch(m_volume_type)
        {
        default:
            {
                mi::Uint32 nb_values = m_brick_size.x*m_brick_size.y*m_brick_size.z;
                for(mi::Uint32 i=0; i<nb_values; i++)
                    brick[i] = (mi::Uint8) (brick_idx*4)%255;
            }
            break;
        
        case TURBULENCE:
            {
                mi::Uint32 brick_idx = 0;
                for(mi::Sint32 x = bbox.min.x; x < bbox.max.x; x++)
                    for(mi::Sint32 y = bbox.min.y; y < bbox.max.y; y++)
                        for(mi::Sint32 z = bbox.min.z; z < bbox.max.z; z++, brick_idx++)
                        {
                            mi::Float32 xf = (mi::Float32) x / (mi::Float32) m_volume_size.x;
                            mi::Float32 yf = (mi::Float32) y / (mi::Float32) m_volume_size.y;
                            mi::Float32 zf = (mi::Float32) z / (mi::Float32) m_volume_size.z;
                            
                            mi::Float32 t = turbulence(20.f*xf, 20.f*yf, 20.f*zf, 1.f);
                            
                            // INFO_LOG << "VALUES: " << bbox << ", " << m_volume_size << ", " << x << ", " << y << ", " << z << ", " << xf << ", " << yf << ", " << zf << ", " << t;
                            brick[brick_idx] = (mi::Uint8) (t*255.f);
                            // INFO_LOG << "brick value: " <<  (int) brick[brick_idx];
                        }
            }
            break;
        }
    }
}


void Distributed_volume_application::get_brick(
    mi::Uint32                      brick_idx, 
    mi::math::Bbox<mi::Float32, 3>   &bbox, 
    mi::Uint8                       *brick) const
{
    mi::Uint32 idx = brick_idx;
    // convert linear brick index to x,y,z position
    // get z
    mi::Uint32 z = idx / (m_volume_bricks_size.x*m_volume_bricks_size.y);
    idx = idx - z*(m_volume_bricks_size.x*m_volume_bricks_size.y);
    
    // get y
    mi::Uint32 y = idx / m_volume_bricks_size.x;

    // get x
    mi::Uint32 x = idx - y*m_volume_bricks_size.x;
    
    //calculate brick bounding box
    bbox.min.x = x*m_brick_size.x;
    bbox.min.y = y*m_brick_size.y;
    bbox.min.z = z*m_brick_size.z;
    
    bbox.max.x = bbox.min.x + m_brick_size.x;
    bbox.max.y = bbox.min.y + m_brick_size.y;
    bbox.max.z = bbox.min.z + m_brick_size.z;
    
    // INFO_LOG << "Brick (" << brick_idx << ") coordinates: " << x << ", " << y << ", " << z <<
        // " and BBox: " << bbox;
        
    // calculate data
    if(brick != NULL)
    {
        switch(m_volume_type)
        {
        default:
            {
                mi::Uint32 nb_values = m_brick_size.x*m_brick_size.y*m_brick_size.z;
                for(mi::Uint32 i=0; i<nb_values; i++)
                    brick[i] = (mi::Uint8) (brick_idx*4)%255;
            }
            break;
        
        case TURBULENCE:
            {
                mi::Uint32 brick_idx = 0;
                for(mi::Sint32 x = bbox.min.x; x < bbox.max.x; x++)
                    for(mi::Sint32 y = bbox.min.y; y < bbox.max.y; y++)
                        for(mi::Sint32 z = bbox.min.z; z < bbox.max.z; z++, brick_idx++)
                        {
                            mi::Float32 xf = (mi::Float32) x / (mi::Float32) m_volume_size.x;
                            mi::Float32 yf = (mi::Float32) y / (mi::Float32) m_volume_size.y;
                            mi::Float32 zf = (mi::Float32) z / (mi::Float32) m_volume_size.z;
                            
                            mi::Float32 t = turbulence(20.f*xf, 20.f*yf, 20.f*zf, 1.f);
                            
                            // INFO_LOG << "VALUES: " << bbox << ", " << m_volume_size << ", " << x << ", " << y << ", " << z << ", " << xf << ", " << yf << ", " << zf << ", " << t;
                            brick[brick_idx] = (mi::Uint8) (t*255.f);
                            // INFO_LOG << "brick value: " <<  (int) brick[brick_idx];
                        }
            }
            break;
        }
    }
}

void Distributed_volume_application::generate_noise()
{
    for(mi::Uint32 x = 0; x < NOISE_SIZE; x++)
        for(mi::Uint32 y = 0; y < NOISE_SIZE; y++)
            for(mi::Uint32 z = 0; z < NOISE_SIZE; z++)
                m_noise[x][y][z] = (rand() % 32768) / 32768.f;
}

mi::Float32 Distributed_volume_application::smooth_noise(
    mi::Float32 x, 
    mi::Float32 y, 
    mi::Float32 z) const
{
   //get fractional part of x and y
   mi::Float32 fractX = x - mi::Sint32(x);
   mi::Float32 fractY = y - mi::Sint32(y);
   mi::Float32 fractZ = z - mi::Sint32(z);   
  
   //wrap around
   mi::Sint32 x1 = (mi::Sint32(x) + NOISE_SIZE) % NOISE_SIZE;
   mi::Sint32 y1 = (mi::Sint32(y) + NOISE_SIZE) % NOISE_SIZE;
   mi::Sint32 z1 = (mi::Sint32(z) + NOISE_SIZE) % NOISE_SIZE;
  
   //neighbor values
   mi::Sint32 x2 = (x1 + NOISE_SIZE - 1) % NOISE_SIZE;
   mi::Sint32 y2 = (y1 + NOISE_SIZE - 1) % NOISE_SIZE;
   mi::Sint32 z2 = (z1 + NOISE_SIZE - 1) % NOISE_SIZE;

   //smooth the noise with bilinear interpolation
   mi::Float32 value = 0.f;
   value += fractX       * fractY       * fractZ       * m_noise[x1][y1][z1];
   value += fractX       * (1.f - fractY) * fractZ       * m_noise[x1][y2][z1];
   value += (1.f - fractX) * fractY       * fractZ       * m_noise[x2][y1][z1];
   value += (1.f - fractX) * (1.f - fractY) * fractZ       * m_noise[x2][y2][z1];

   value += fractX       * fractY       * (1.f - fractZ) * m_noise[x1][y1][z2];
   value += fractX       * (1.f - fractY) * (1.f - fractZ) * m_noise[x1][y2][z2];
   value += (1.f - fractX) * fractY       * (1.f - fractZ) * m_noise[x2][y1][z2];
   value += (1.f - fractX) * (1.f - fractY) * (1.f - fractZ) * m_noise[x2][y2][z2];

   return value;    
}
    
mi::Float32 Distributed_volume_application::turbulence(
    mi::Float32 x, 
    mi::Float32 y, 
    mi::Float32 z, 
    mi::Float32 size) const
{
mi::Float32 value = 0.f, initial_size = size;
   
    while(size >= 1.f)
    {
        value += smooth_noise(x / size, y / size, z / size) * size;
        size /= 2.f;
    }
   
    return(value / initial_size);    
}
