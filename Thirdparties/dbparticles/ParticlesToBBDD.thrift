/*
 * BigDataParicles interface
 */

namespace cpp BBDDParticles
namespace java es.bsc.sparkandra.BBDDparticles

//Tipo para las respuestas, output
struct Particle {
  1: double time,
  2: i32 partId,
  3: double x,
  4: double y,
  5: double z,
  6: map<string,double> doubleProperties
  7: map<string,i32> intProperties
}

struct InsertData{
    1:list<Particle> data,
    2:i32 time
}

struct InsertResult{
    1:bool sucess,
    2:string error
}

service ParticlesToBBDD {

  // Insert statement
  InsertResult insertParticles(1:InsertData data),

  //solo para test
  i32 ping(),

}
