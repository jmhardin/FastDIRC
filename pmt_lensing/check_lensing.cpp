#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <math.h>

#include <Goptical/Sys/System>
#include <Goptical/Analysis/Spot>

#include <Goptical/Math/Vector>
#include <Goptical/Math/VectorPair>

#include <Goptical/Sys/Element>
#include <Goptical/Sys/Container>
#include <Goptical/Sys/System>
#include <Goptical/Sys/SourcePoint>
#include <Goptical/Sys/SourceRays>
#include <Goptical/Sys/Source>
#include <Goptical/Sys/Group>
#include <Goptical/Design/Telescope/Newton>
#include <Goptical/Sys/Image>
#include <Goptical/Sys/Lens>
#include <Goptical/Sys/OpticalSurface>
#include <Goptical/Sys/Mirror>
#include <Goptical/Sys/Mirror>
#include <Goptical/Data/Plot>
#include <Goptical/Data/PlotData>
#include <Goptical/Data/Set>

#include <Goptical/Trace/Tracer>
#include <Goptical/Trace/Result>
#include <Goptical/Trace/Distribution>
#include <Goptical/Trace/Sequence>
#include <Goptical/Trace/Params>
#include <Goptical/Trace/Ray>

#include <Goptical/Light/Ray>
#include <Goptical/Light/SpectralLine>
#include "Goptical/common.hh"

#include <Goptical/Io/RendererSvg>
#include <Goptical/Io/Rgb>
#include <Goptical/Io/RendererViewport>

#include <Goptical/Material/SellmeierMod>

#include <Goptical/Math/Transform>

#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TRandom3.h>

#include <stdlib.h>
#include <math.h>

#include <TRandom3.h>

#include "../include/dirc_optical_components.h"

//This is just a toy sim, so I'm keeping it all in one file
struct vec3
{
  double x;
  double y;
  double z;
};

void gen_random_ray(\
  double source_zoff,\
  double source_half_height,\
  double lens_half_height,\
  TRandom3 *rand_gen,\
  vec3 &rpos,\
  vec3 &rdir)
{
  double start_y = rand_gen->Uniform(-source_half_height,source_half_height);
  double target_y = rand_gen->Uniform(-lens_half_height,lens_half_height);
  
  rpos.x = 0;
  rpos.y = start_y;
  rpos.z = source_zoff;

  double y_diff = target_y - start_y;
  double dir_r = sqrt(y_diff*y_diff + source_zoff*source_zoff);
  
  rdir.x = 0;
  rdir.y = y_diff/dir_r;
  rdir.z = -source_zoff/dir_r;
  
}
  
void trace_rays(\
  Sys::System &sys,\
  Sys::SourceRays &srcrays, \
  Sys::Image &image,\
  bool outframe,\
  vec3 &rval)
{
  
  sys.add(srcrays);
  
  Trace::Tracer tracer(sys);
  tracer.get_params().set_max_bounce(5000);
  tracer.get_trace_result().set_generated_save_state(srcrays);
  tracer.get_trace_result().set_intercepted_save_state(image);

  tracer.trace();

  _Goptical::Trace::rays_queue_t rays = tracer.get_trace_result().get_intercepted(image);

  double x, y;
  
  for (unsigned int i = 0; i < rays.size(); i++)
  {
	  x = ((Math::Vector3)rays[i]->get_intercept_point()).x();
	  y = ((Math::Vector3)rays[i]->get_intercept_point()).y();
	  
	  rval.x = x;
	  rval.y = y;
	  
  }
//   printf("%d\n",rays.size());

  if (outframe == true)
  {
	  printf("Rendering 3d\n");
	  
	  Io::RendererSvg       svg_renderer("layout.svg", 1000, 700);
	  Io::RendererViewport  &renderer = svg_renderer;
	  // 3d system layout
	  
// 	  renderer.set_feature_size(200);
	  
	  renderer.set_perspective();
// 	  renderer.set_fov(3);	
	  
	  
	  sys.draw_3d_fit(renderer, 0);
// 	  renderer.set_camera_transform(renderer.get_camera_transform().linear_rotation(Math::Vector3(0,0,0)));
	  sys.draw_3d(renderer);
	  
	  tracer.get_trace_result().draw_3d(renderer);
	  
	  printf("rendered 3d\n");
  }

  sys.remove(srcrays);

}


vec3 propagate_light(\
  double r,\
  double lens_height,\
  double zoff_target,\
  double isystem_width,\
  vec3 ipos,\
  vec3 idir,\
  int num_display = 0,\
  double zoff_display = 0,\
  double yheight_display = 0)
{
  vec3 rval;
  
  double nl = 1.33;
  double nq = 1.47;
  
  double system_width = isystem_width;
  bool out_layout = false;
  
  Sys::System sys;
  
//   ref<qwartzCrystal> qwartz = ref<qwartzCrystal>::create();
//   ref<boxLiquid> liquid = ref<boxLiquid>::create();
  
  Material::Sellmeier liquid(0,0,0,0,0,0);
  Material::Sellmeier qwartz(0,0,0,0,0,0);
  
  //Actually mispelled in gopticals header :(
  liquid.set_contant_term(nl*nl);
  qwartz.set_contant_term(nq*nq);
  
  Sys::SourceRays srcrays(Math::Vector<3>(0,0,0));
  
  srcrays.set_material(liquid);
//   srcrays.add_spectral_line(Light::SpectralLine(550));
  
  
  if (num_display == 0)
  {
    Math::VectorPair3 new_ray(\
	  ipos.x,\
	  ipos.y,\
	  ipos.z,\
	  idir.x,\
	  idir.y,\
	  idir.z);
    
    srcrays.add_rays(new_ray,&srcrays);
  }
  else
  {
    TRandom3 *rgen = new TRandom3();
    for (int i = 0; i < num_display; i++)
    {
      vec3 pos;
      vec3 dir;
      
      gen_random_ray(\
	-zoff_display,\
	yheight_display/2,\
	lens_height/2,\
	rgen,\
	pos,\
	dir);
      
      Math::VectorPair3 new_ray(\
	pos.x,\
	pos.y,\
	pos.z,\
	dir.x,\
	dir.y,\
	dir.z);
    
      
//       printf("%12.04f %12.04f %12.04f %12.04f %12.04f %12.04f \n",pos.x,pos.y,pos.z,dir.x,dir.y,dir.z);
      
      
      srcrays.add_rays(new_ray,&srcrays);
    }
    out_layout = true;
  }
      
    
    
//   sys.add(srcrays);
  
  
  ref<Curve::Flat> image_curve = ref<Curve::Flat>::create();
  ref<Shape::Rectangle> image_shape = ref<Shape::Rectangle>::create(system_width,3*lens_height);
  
  Sys::Image image(\
      Math::Vector<3>(0,0,zoff_target),\
      image_curve,\
      image_shape);
  
  sys.add(image);
  
  if (out_layout == true)
  {
    Sys::Mirror display_help_image(\
      Math::Vector<3>(0,0,-1.1*zoff_display),\
      image_curve,\
      image_shape);
    
//     sys.add(display_help_image);
  }

  
  ref<cylinderY> lens_curve = ref<cylinderY>::create(r);
  ref<Shape::Rectangle> lens_shape = ref<Shape::Rectangle>::create(system_width, lens_height);
  
  double sagitta_off = r - sqrt(r*r - lens_height*lens_height/4);
  
  ref<Sys::OpticalSurface> lens_surface = \
    ref<Sys::OpticalSurface>::create(\
      Math::VectorPair<3>(0,0,-sagitta_off),\
      lens_curve,\
      lens_shape,\
      qwartz,\
      liquid);
  
  lens_surface->rotate(0,180,0);
    
  sys.add(*lens_surface);
  
//   Sys::OpticalSurface lens_surface(\
//     Math::VectorPair<3>(0,0,-sagitta_off),\
//     -r,\
//     lens_height/2,\
//     liquid,\
//     qwartz);
//   
//   sys.add(lens_surface);
  
  trace_rays(\
    sys,\
    srcrays, \
    image,\
    out_layout,\
    rval);
  
  rval.z = zoff_target;
  
  return rval;
}
int main(int argc, char** argv)
{
  int num_points_tested = 1000000;
  int num_points_displayed = 100;
  
  
  double pmt_height = 6;
  double lens_r = 50;
  double lens_height = 10;
  double zoff_target = 0;
  double zoff_source = 5*lens_height;
  double system_width = lens_height;
  double source_height = 2*zoff_source;
  
  vec3 pos,dir;
  
  propagate_light(\
    lens_r,\
    lens_height,\
    zoff_target,\
    system_width,\
    pos,\
    dir,\
    num_points_displayed,\
    zoff_source,\
    source_height);
  
  TRandom3 *rgen = new TRandom3;
  
  char* rootfilename = new char[256];
  sprintf(rootfilename,"pmt_lens.root");
  
  TFile* tfile = new TFile(rootfilename,"RECREATE");
  
  TH2F* produced_rays = new TH2F("produced_rays","Location/direction of Produced Rays",1000,-source_height/2,source_height/2,1000,-3.14,3.14);
  TH2F* detected_rays = new TH2F("detected_rays","Initial Location/direction of Detected Rays",1000,-source_height/2,source_height/2,1000,-3.14,3.14);
  
  double theta, starty;
  vec3 rpos;
  int detected = 0;
  
  for (int i = 0; i < num_points_tested; i++)
  {
    
    gen_random_ray(\
      zoff_source,\
      source_height/2,\
      lens_height/2,\
      rgen,\
      pos,\
      dir);
    
    theta = atan(dir.y/dir.z);
    starty = pos.y;
    
    produced_rays->Fill(starty,theta);
    
    dir.z = 1;
    dir.y = 0;
    pos.y = -1;
    
    rpos = propagate_light(\
      lens_r,\
      lens_height,\
      zoff_target,\
      system_width,\
      pos,\
      dir);
    
    printf("%12.04f\n",rpos.y);
    
    if (fabs(rpos.y) < pmt_height/2)
    {
      detected_rays->Fill(starty,theta);
      detected++;
    }
    
  }
  
  printf("Fraction detected: %d/%d\n",detected,num_points_tested);
  produced_rays->Write();
  detected_rays->Write();
  
  tfile->Close();
  
  return 0;
}