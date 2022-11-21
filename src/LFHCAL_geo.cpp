//==========================================================================
//  Implementation of longitudinally separated forward calorimeter
//--------------------------------------------------------------------------
//  Author: Friederike Bock (ORNL)
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>
#include "XML/Layering.h"
#include "XML/Utilities.h"
using namespace dd4hep;

struct moduleParamsStrct{
  moduleParamsStrct(): mod_BIwidth(0.), mod_BIheight(0.), mod_SWThick(0.),   mod_TWThick(0.), mod_FWThick (0.),
                      mod_BWThick(0.), mod_width(0.), mod_height(0.), mod_nodgeWidthAbsA(0.), mod_nodgeWidthAbsB(0.), mod_nodgeWidthAbsC (0.),
                      mod_nodgeWidthScin(0.), mod_nodgeDepth(0.), mod_sepDepth(0.),  mod_visStr(""), mod_regStr(""), mod_limStr("")
                      {}
  moduleParamsStrct(   double BIwidth, double BIheight, double SWThick, double TWThick, double FWThick, double BWThick, double width, double height, 
                       double nodgeWidthAbsA, double nodgeWidthAbsB, double nodgeWidthAbsC, double nodgeWidthScin, double nodgeDepth, double sepDepth, 
                       std::string visStr, std::string regStr, std::string limStr){
      mod_BIwidth       = BIwidth;
      mod_BIheight      = BIheight;
      mod_SWThick       = SWThick;  
      mod_TWThick       = TWThick;
      mod_FWThick       = FWThick;
      mod_BWThick       = BWThick;
      mod_width         = width;
      mod_height        = height;
      mod_nodgeWidthAbsA  = nodgeWidthAbsA;
      mod_nodgeWidthAbsB  = nodgeWidthAbsB;
      mod_nodgeWidthAbsC  = nodgeWidthAbsC;
      mod_nodgeWidthScin  = nodgeWidthScin;
      mod_nodgeDepth      = nodgeDepth;
      mod_sepDepth        = sepDepth;
      mod_visStr          = visStr;
      mod_regStr          = regStr;
      mod_limStr          = limStr;
  }
  double      mod_BIwidth;
  double      mod_BIheight;
  double      mod_SWThick;
  double      mod_TWThick;
  double      mod_FWThick;
  double      mod_BWThick;
  double      mod_width;
  double      mod_height;
  double      mod_nodgeWidthAbsA;
  double      mod_nodgeWidthAbsB;
  double      mod_nodgeWidthAbsC;
  double      mod_nodgeWidthScin;
  double      mod_nodgeDepth;
  double      mod_sepDepth;
  std::string mod_visStr;
  std::string mod_regStr;
  std::string mod_limStr;
} ;

struct sliceParamsStrct{
  sliceParamsStrct(): layer_ID(0), slice_ID(0), slice_partID(0), slice_thick(0.), slice_offset(0.), slice_matStr(""), slice_visStr(""), slice_regStr(""), slice_limStr("")
                      {}
  sliceParamsStrct(                   
                      int l_ID, int sl_ID, int sl_partID, double sl_thick, double sl_off, std::string sl_matStr, std::string sl_visStr, std::string sl_regStr, std::string sl_limStr ){
      layer_ID      = l_ID;
      slice_ID      = sl_ID;
      slice_partID  = sl_partID;
      slice_thick   = sl_thick;
      slice_offset  = sl_off;
      slice_matStr  = sl_matStr;
      slice_visStr  = sl_visStr;
      slice_regStr  = sl_regStr;
      slice_limStr  = sl_limStr;
      
  }
  int         layer_ID;
  int         slice_ID;
  int         slice_partID;
  double      slice_thick;
  double      slice_offset;
  std::string slice_matStr;
  std::string slice_visStr;
  std::string slice_regStr;
  std::string slice_limStr;
};  


//************************************************************************************************************
//************************** Assembly for absorber plates for 8M modules *************************************
//************************************************************************************************************
Assembly createAbsorberPlateEightM(Detector& desc,
                                   std::string basename, 
                                   double h_mod,
                                   double w_mod,
                                   double t_mod_tp,
                                   double t_mod_sp,
                                   double t_slice,
                                   double h_nodge,
                                   double w_nodgeA,
                                   double w_nodgeB,
                                   double w_nodgeC,
                                   Material slice_mat,
                                   std::string region,
                                   std::string limit,
                                   std::string vis
){
  Assembly modAbsAssembly(basename);
          
  Box         slice_plate((w_mod-2*t_mod_sp) / 2., (h_mod-2*t_mod_tp-2*h_nodge) / 2., t_slice / 2.);
  Box         nodge_1(w_nodgeA / 2., h_nodge / 2., t_slice / 2.);
  Box         nodge_2(w_nodgeB / 2., h_nodge / 2., t_slice / 2.);
  Box         nodge_3(w_nodgeC / 2., h_nodge / 2., t_slice / 2.);
  
  Volume      slice_vol(basename+"main", slice_plate, slice_mat);
  Volume      nodge_t1_vol(basename+"_nodge_t1", nodge_3, slice_mat);
  Volume      nodge_t2_vol(basename+"_nodge_t2", nodge_1, slice_mat);
  Volume      nodge_t3_vol(basename+"_nodge_t3", nodge_2, slice_mat);
  Volume      nodge_t4_vol(basename+"_nodge_t4", nodge_1, slice_mat);
  Volume      nodge_t5_vol(basename+"_nodge_t5", nodge_3, slice_mat);
  Volume      nodge_b1_vol(basename+"_nodge_b1", nodge_3, slice_mat);
  Volume      nodge_b2_vol(basename+"_nodge_b2", nodge_1, slice_mat);
  Volume      nodge_b3_vol(basename+"_nodge_t3", nodge_2, slice_mat);
  Volume      nodge_b4_vol(basename+"_nodge_b4", nodge_1, slice_mat);
  Volume      nodge_b5_vol(basename+"_nodge_b5", nodge_3, slice_mat);
  // Setting slice attributes
  slice_vol.setAttributes(desc, region, limit, vis);
  nodge_t1_vol.setAttributes(desc, region, limit, vis);
  nodge_t2_vol.setAttributes(desc, region, limit, vis);
  nodge_t3_vol.setAttributes(desc, region, limit, vis);
  nodge_t4_vol.setAttributes(desc, region, limit, vis);
  nodge_t5_vol.setAttributes(desc, region, limit, vis);
  nodge_b1_vol.setAttributes(desc, region, limit, vis);
  nodge_b2_vol.setAttributes(desc, region, limit, vis);
  nodge_b3_vol.setAttributes(desc, region, limit, vis);
  nodge_b4_vol.setAttributes(desc, region, limit, vis);
  nodge_b5_vol.setAttributes(desc, region, limit, vis);
  
  // Placing slice within slice
  PlacedVolume pvm;
  // base plate
  pvm = modAbsAssembly.placeVolume(slice_vol, Position(0, 0, 0 ));

  // top row of nodges
  pvm = modAbsAssembly.placeVolume(nodge_t1_vol, Position(w_mod/2-t_mod_sp-0.5*w_nodgeC, h_mod/2-t_mod_tp-0.5*h_nodge, 0 ));  
  pvm = modAbsAssembly.placeVolume(nodge_t2_vol, Position((w_mod-2*t_mod_sp)/4, h_mod/2-t_mod_tp-0.5*h_nodge, 0 ));
  pvm = modAbsAssembly.placeVolume(nodge_t3_vol, Position(0, h_mod/2-t_mod_tp-0.5*h_nodge, 0 ));
  pvm = modAbsAssembly.placeVolume(nodge_t4_vol, Position(-(w_mod-2*t_mod_sp)/4, h_mod/2-t_mod_tp-0.5*h_nodge, 0 ));
  pvm = modAbsAssembly.placeVolume(nodge_t5_vol, Position(-w_mod/2+t_mod_sp+0.5*w_nodgeC, h_mod/2-t_mod_tp-0.5*h_nodge, 0 ));
  
  // bottom row of nodges
  pvm = modAbsAssembly.placeVolume(nodge_b1_vol, Position(w_mod/2-t_mod_sp-0.5*w_nodgeC, -(h_mod-2*t_mod_tp)/2+0.5*h_nodge, 0 ));
  pvm = modAbsAssembly.placeVolume(nodge_b2_vol, Position((w_mod-2*t_mod_sp)/4, -(h_mod-2*t_mod_tp)/2+0.5*h_nodge, 0 ));
  pvm = modAbsAssembly.placeVolume(nodge_b3_vol, Position(0, -(h_mod-2*t_mod_tp)/2 +0.5*h_nodge, 0 ));
  pvm = modAbsAssembly.placeVolume(nodge_b4_vol, Position(-(w_mod-2*t_mod_sp)/4, -(h_mod-2*t_mod_tp)/2 +0.5*h_nodge, 0 ));
  pvm = modAbsAssembly.placeVolume(nodge_b5_vol, Position(-w_mod/2+t_mod_sp+0.5*w_nodgeC, -(h_mod-2*t_mod_tp)/2 +0.5*h_nodge, 0 ));
  
  return modAbsAssembly;
}

//************************************************************************************************************
//************************** Assembly for absorber plates for 4M modules *************************************
//************************************************************************************************************
Assembly createAbsorberPlateFourM(Detector& desc,
                                   std::string basename, 
                                   double h_mod,
                                   double w_mod,
                                   double t_mod_tp,
                                   double t_mod_sp,
                                   double t_slice,
                                   double h_nodge,
                                   double w_nodgeA,
                                   double w_nodgeC,
                                   Material slice_mat,
                                   std::string region,
                                   std::string limit,
                                   std::string vis
){
  Assembly modAbsAssembly(basename);
          
  Box         slice_plate((w_mod-2*t_mod_sp) / 2., (h_mod-2*t_mod_tp-2*h_nodge) / 2., t_slice / 2.);
  Box         nodge_1(w_nodgeA / 2., h_nodge / 2., t_slice / 2.);
  Box         nodge_3(w_nodgeC / 2., h_nodge / 2., t_slice / 2.);
  
  Volume      slice_vol(basename+"main", slice_plate, slice_mat);
  Volume      nodge_t1_vol(basename+"_nodge_t1", nodge_3, slice_mat);
  Volume      nodge_t2_vol(basename+"_nodge_t2", nodge_1, slice_mat);
  Volume      nodge_t3_vol(basename+"_nodge_t3", nodge_3, slice_mat);
  Volume      nodge_b1_vol(basename+"_nodge_b1", nodge_3, slice_mat);
  Volume      nodge_b2_vol(basename+"_nodge_b2", nodge_1, slice_mat);
  Volume      nodge_b3_vol(basename+"_nodge_b3", nodge_3, slice_mat);
  // Setting slice attributes
  slice_vol.setAttributes(desc, region, limit, vis);
  nodge_t1_vol.setAttributes(desc, region, limit, vis);
  nodge_t2_vol.setAttributes(desc, region, limit, vis);
  nodge_t3_vol.setAttributes(desc, region, limit, vis);
  nodge_b1_vol.setAttributes(desc, region, limit, vis);
  nodge_b2_vol.setAttributes(desc, region, limit, vis);
  nodge_b3_vol.setAttributes(desc, region, limit, vis);
  
  // Placing slice within slice
  PlacedVolume pvm;
  // base plate
  pvm = modAbsAssembly.placeVolume(slice_vol, Position(0, 0, 0 ));

  // top row of nodges
  pvm = modAbsAssembly.placeVolume(nodge_t1_vol, Position(w_mod/2-t_mod_sp-0.5*w_nodgeC, h_mod/2-t_mod_tp-0.5*h_nodge, 0 ));  
  pvm = modAbsAssembly.placeVolume(nodge_t2_vol, Position(0, h_mod/2-t_mod_tp-0.5*h_nodge, 0 ));
  pvm = modAbsAssembly.placeVolume(nodge_t3_vol, Position(-w_mod/2+t_mod_sp+0.5*w_nodgeC, h_mod/2-t_mod_tp-0.5*h_nodge, 0 ));
  
  // bottom row of nodges
  pvm = modAbsAssembly.placeVolume(nodge_b1_vol, Position(w_mod/2-t_mod_sp-0.5*w_nodgeC, -(h_mod-2*t_mod_tp)/2+0.5*h_nodge, 0 ));
  pvm = modAbsAssembly.placeVolume(nodge_b2_vol, Position(0, -(h_mod-2*t_mod_tp)/2 +0.5*h_nodge, 0 ));
  pvm = modAbsAssembly.placeVolume(nodge_b3_vol, Position(-w_mod/2+t_mod_sp+0.5*w_nodgeC, -(h_mod-2*t_mod_tp)/2 +0.5*h_nodge, 0 ));
  
  return modAbsAssembly;
}


//************************************************************************************************************
//************************** Filler plate i.e. Tyvek/Air for 8M module ***************************************
//************************************************************************************************************
Assembly createFillerPlateEightM( Detector& desc,
                                   std::string basename, 
                                   double h_mod,
                                   double w_mod,
                                   double t_mod_tp,
                                   double t_mod_sp,
                                   double t_slice,
                                   double h_nodge,
                                   double w_nodge,
                                   Material slice_mat,
                                   std::string region,
                                   std::string limit,
                                   std::string vis
){
  Assembly modFillAssembly(basename);
  
  Box         slice_plate((w_mod-2*t_mod_sp) / 2., (h_mod-2*t_mod_tp-2*h_nodge) / 2., t_slice / 2.);
  Box         nodge_1(w_nodge / 2., h_nodge / 2., t_slice / 2.);
  
  Volume      slice_vol(basename+"main", slice_plate, slice_mat);
  Volume      nodge_t1_vol(basename+"_nodge_t1", nodge_1, slice_mat);
  Volume      nodge_t2_vol(basename+"_nodge_t2", nodge_1, slice_mat);
  Volume      nodge_b1_vol(basename+"_nodge_b1", nodge_1, slice_mat);
  Volume      nodge_b2_vol(basename+"_nodge_b2", nodge_1, slice_mat);
  
  // Setting slice attributes
  slice_vol.setAttributes(desc, region, limit, vis);
  nodge_t1_vol.setAttributes(desc, region, limit, vis);
  nodge_t2_vol.setAttributes(desc, region, limit, vis);
  nodge_b1_vol.setAttributes(desc, region, limit, vis);
  nodge_b2_vol.setAttributes(desc, region, limit, vis);
  
  // Placing slice within slice
  PlacedVolume pvm;
  pvm = modFillAssembly.placeVolume(slice_vol, Position(0, 0, 0 ));
  pvm = modFillAssembly.placeVolume(nodge_t1_vol, Position((w_mod-2*t_mod_sp)/4, h_mod/2-t_mod_tp-0.5*h_nodge, 0 ));
  pvm = modFillAssembly.placeVolume(nodge_t2_vol, Position(-(w_mod-2*t_mod_sp)/4, h_mod/2-t_mod_tp-0.5*h_nodge, 0 ));
  pvm = modFillAssembly.placeVolume(nodge_b1_vol, Position((w_mod-2*t_mod_sp)/4, -(h_mod-2*t_mod_tp)/2+0.5*h_nodge, 0 ));
  pvm = modFillAssembly.placeVolume(nodge_b2_vol, Position(-(w_mod-2*t_mod_sp)/4, -(h_mod-2*t_mod_tp)/2 +0.5*h_nodge, 0 ));
  
  return modFillAssembly;
}

//************************************************************************************************************
//************************** Filler plate i.e. Tyvek/Air for 4M module ***************************************
//************************************************************************************************************
Assembly createFillerPlateFourM( Detector& desc,
                                   std::string basename, 
                                   double h_mod,
                                   double w_mod,
                                   double t_mod_tp,
                                   double t_mod_sp,
                                   double t_slice,
                                   double h_nodge,
                                   double w_nodge,
                                   Material slice_mat,
                                   std::string region,
                                   std::string limit,
                                   std::string vis
){
  Assembly modFillAssembly(basename);
  
  Box         slice_plate((w_mod-2*t_mod_sp) / 2., (h_mod-2*t_mod_tp-2*h_nodge) / 2., t_slice / 2.);
  Box         nodge_1(w_nodge / 2., h_nodge / 2., t_slice / 2.);
  
  Volume      slice_vol(basename+"main", slice_plate, slice_mat);
  Volume      nodge_t1_vol(basename+"_nodge_t1", nodge_1, slice_mat);
  Volume      nodge_b1_vol(basename+"_nodge_b1", nodge_1, slice_mat);
  
  // Setting slice attributes
  slice_vol.setAttributes(desc, region, limit, vis);
  nodge_t1_vol.setAttributes(desc, region, limit, vis);
  nodge_b1_vol.setAttributes(desc, region, limit, vis);
  
  // Placing slice within slice
  PlacedVolume pvm;
  pvm = modFillAssembly.placeVolume(slice_vol, Position(0, 0, 0 ));
  pvm = modFillAssembly.placeVolume(nodge_t1_vol, Position(0, h_mod/2-t_mod_tp-0.5*h_nodge, 0 ));
  pvm = modFillAssembly.placeVolume(nodge_b1_vol, Position(0, -(h_mod-2*t_mod_tp)/2+0.5*h_nodge, 0 ));
  
  return modFillAssembly;
}


//************************************************************************************************************
//************************** single scintillator plate for tower *********************************************
//************************************************************************************************************
Assembly createScintillatorTower( Detector& desc,
                                   std::string basename, 
                                   double w_tow,
                                   double h_tow,
                                   double t_slice,
                                   double h_nodge,
                                   double w_nodge,
                                   Material slice_mat,
                                   std::string region,
                                   std::string limit,
                                   std::string vis,
                                  SensitiveDetector sens
){
  
  Assembly towerAss(basename);
          
  Box         slice_plate(w_tow / 2., h_tow / 2., t_slice / 2.);
  Box         slice_nodge( w_nodge/ 4., h_nodge / 2., t_slice / 2.);

  Volume      slice_vol(basename+"_main", slice_plate, slice_mat);
  Volume      slicenodge_vol(basename+"nodge", slice_nodge, slice_mat);
  
  // Setting appropriate slices as sensitive
  sens.setType("calorimeter");
  slice_vol.setSensitiveDetector(sens);
  slicenodge_vol.setSensitiveDetector(sens);

  // Setting slice attributes
  slice_vol.setAttributes(desc, region, limit, vis);
  slicenodge_vol.setAttributes(desc, region, limit, vis);
  // Placing slice within tower slice
  PlacedVolume pvm;
  pvm = towerAss.placeVolume(slice_vol, Position(0., 0., 0));
  pvm.addPhysVolID("part",0);
  pvm = towerAss.placeVolume(slicenodge_vol, Position(w_tow/2-w_nodge/ 4., h_tow/2+0.5*h_nodge, 0));
  pvm.addPhysVolID("part",1);
  return towerAss;
}  

//************************************************************************************************************
//************************** create scintillator plate with separations for 8M *******************************
//************************************************************************************************************
Assembly createScintillatorPlateEightM( Detector& desc,
                                        std::string basename, 
                                        int modID,
                                        int layerID,
                                        double h_mod,
                                        double w_mod,
                                        double t_mod_tp,
                                        double t_mod_sp,
                                        double t_slice,
                                        double h_nodge,
                                        double w_nodge,
                                        double t_sep,
                                        Material slice_mat,
                                        std::string region,
                                        std::string limit,
                                        std::string vis,
                                        SensitiveDetector sens
){
  // Tower placement in 8M module
  //======================================================================
  //||              ||              ||                ||                ||
  //||      0       ||      1       ||        2       ||        3       ||
  //||              ||              ||                ||                ||
  //======================================================================
  //||              ||              ||                ||                ||
  //||      4       ||      5       ||        6       ||        7       ||
  //||              ||              ||                ||                ||
  //======================================================================
  Assembly modScintAssembly(basename);
  double towerWidth   = (w_mod-2*t_mod_sp-3*t_sep)/4;
  double towerHeight  = (h_mod-2*t_mod_tp-t_sep-2*h_nodge)/2;

  // placement volumes
  PlacedVolume pvm;
  
  // titanium-dioxide separations
  Box         slice_ti02hor( (w_mod-2*t_mod_sp) / 2., t_sep / 2., t_slice / 2.);
  Box         slice_ti02ver1( t_sep / 2., (towerHeight+h_nodge)/2 , t_slice / 2.);
  Box         slice_ti02ver2( t_sep / 2., (towerHeight)/2 , t_slice / 2.);
  Volume      sliceti02H_vol(basename+"ti02_h", slice_ti02hor, desc.material("Ti02Epoxy"));
  Volume      sliceti02t1_vol(basename+"ti02_t1", slice_ti02ver1, desc.material("Ti02Epoxy"));
  Volume      sliceti02t2_vol(basename+"ti02_t2", slice_ti02ver2, desc.material("Ti02Epoxy"));
  Volume      sliceti02t3_vol(basename+"ti02_t3", slice_ti02ver1, desc.material("Ti02Epoxy"));
  Volume      sliceti02b1_vol(basename+"ti02_b1", slice_ti02ver1, desc.material("Ti02Epoxy"));
  Volume      sliceti02b2_vol(basename+"ti02_b2", slice_ti02ver2, desc.material("Ti02Epoxy"));
  Volume      sliceti02b3_vol(basename+"ti02_b3", slice_ti02ver1, desc.material("Ti02Epoxy"));

  sliceti02H_vol.setAttributes(desc, region, limit, "LFHCALLayerTiOVis");
  sliceti02t1_vol.setAttributes(desc, region, limit,"LFHCALLayerTiOVis");
  sliceti02t2_vol.setAttributes(desc, region, limit, "LFHCALLayerTiOVis");
  sliceti02t3_vol.setAttributes(desc, region, limit, "LFHCALLayerTiOVis");
  sliceti02b1_vol.setAttributes(desc, region, limit, "LFHCALLayerTiOVis");
  sliceti02b2_vol.setAttributes(desc, region, limit, "LFHCALLayerTiOVis");
  sliceti02b3_vol.setAttributes(desc, region, limit, "LFHCALLayerTiOVis");
  
  pvm = modScintAssembly.placeVolume(sliceti02H_vol, Position(0, 0, 0 ));
  pvm = modScintAssembly.placeVolume(sliceti02t1_vol, Position(towerWidth+t_sep, 0.5*(towerHeight+t_sep+h_nodge), 0 ));
  pvm = modScintAssembly.placeVolume(sliceti02t2_vol, Position(0, 0.5*(towerHeight+t_sep), 0 ));
  pvm = modScintAssembly.placeVolume(sliceti02t3_vol, Position(-(towerWidth+t_sep), 0.5*(towerHeight+t_sep+h_nodge), 0 ));
  pvm = modScintAssembly.placeVolume(sliceti02b1_vol, Position(towerWidth+t_sep, -0.5*(towerHeight+t_sep+h_nodge), 0 ));
  pvm = modScintAssembly.placeVolume(sliceti02b2_vol, Position(0, -0.5*(towerHeight+t_sep), 0 ));
  pvm = modScintAssembly.placeVolume(sliceti02b3_vol, Position(-(towerWidth+t_sep), -0.5*(towerHeight+t_sep+h_nodge), 0 ));

  // 8M module placement of scintillator for tower
  double rotZ[8] = {0,    0, 0,     0,  0,    0,    0,    0};
  double rotY[8] = {M_PI, 0, M_PI,  0,  M_PI, 0,    M_PI, 0};
  double rotX[8] = {0,    0, 0,     0,  M_PI, M_PI, M_PI, M_PI};
  double posX[8] = {(towerWidth*1.5+1.5*t_sep),     (towerWidth*0.5+0.5*t_sep),     -(towerWidth*0.5+0.5*t_sep),    -(towerWidth*1.5+1.5*t_sep),
                    (towerWidth*1.5+1.5*t_sep),     (towerWidth*0.5+0.5*t_sep),     -(towerWidth*0.5+0.5*t_sep),    -(towerWidth*1.5+1.5*t_sep)};
  double posY[8] = {0.5*(towerHeight)+0.5*t_sep,    0.5*(towerHeight)+0.5*t_sep,    0.5*(towerHeight)+0.5*t_sep,    0.5*(towerHeight)+0.5*t_sep,
                    -(0.5*(towerHeight)+0.5*t_sep), -(0.5*(towerHeight)+0.5*t_sep), -(0.5*(towerHeight)+0.5*t_sep), -(0.5*(towerHeight)+0.5*t_sep)};
  double posZ[8] = {0,                                        0,                                        0,                                        0,
                    0,                                        0,                                        0,                                        0};
  // loop over all towers within same module
  for (int i = 0; i < 8; i++){
    Assembly modScintTowerAss = createScintillatorTower( desc,  basename+ _toString(i, "_tower_%d"),  
                                                            towerWidth, towerHeight, t_slice,
                                                            h_nodge, w_nodge, 
                                                            slice_mat, region, limit, vis, sens);
    pvm = modScintAssembly.placeVolume(modScintTowerAss, Transform3D(RotationZYX(rotZ[i], rotY[i], rotX[i]), Position(posX[i], posY[i], posZ[i] )));
    pvm.addPhysVolID("module", modID).addPhysVolID("tower", i).addPhysVolID("layer", layerID);
  }
  
  return modScintAssembly;         
}          

//************************************************************************************************************
//************************** create scintillator plate with separations for 4M *******************************
//************************************************************************************************************
Assembly createScintillatorPlateFourM( Detector& desc,
                                        std::string basename, 
                                        int modID,
                                        int layerID,
                                        double h_mod,
                                        double w_mod,
                                        double t_mod_tp,
                                        double t_mod_sp,
                                        double t_slice,
                                        double h_nodge,
                                        double w_nodge,
                                        double t_sep,
                                        Material slice_mat,
                                        std::string region,
                                        std::string limit,
                                        std::string vis,
                                        SensitiveDetector sens
){
  // Tower placement in 4M module
  //==================================
  //||              ||              ||
  //||      0       ||      1       ||
  //||              ||              ||
  //==================================
  //||              ||              ||
  //||      2       ||      3       ||
  //||              ||              ||
  //==================================
  Assembly modScintAssembly(basename);
  double towerWidth   = (w_mod-2*t_mod_sp-t_sep)/2;
  double towerHeight  = (h_mod-2*t_mod_tp-t_sep-2*h_nodge)/2;

  // placement volumes
  PlacedVolume pvm;
  
  // titanium-dioxide separations
  Box         slice_ti02hor( (w_mod-2*t_mod_sp) / 2., t_sep / 2., t_slice / 2.);
  Box         slice_ti02ver1( t_sep / 2., (towerHeight+h_nodge)/2 , t_slice / 2.);
  Volume      sliceti02H_vol(basename+"ti02_h", slice_ti02hor, desc.material("Ti02Epoxy"));
  Volume      sliceti02t_vol(basename+"ti02_t", slice_ti02ver1, desc.material("Ti02Epoxy"));
  Volume      sliceti02b_vol(basename+"ti02_b", slice_ti02ver1, desc.material("Ti02Epoxy"));

  sliceti02H_vol.setAttributes(desc, region, limit, "LFHCALLayerTiOVis");
  sliceti02t_vol.setAttributes(desc, region, limit,"LFHCALLayerTiOVis");
  sliceti02b_vol.setAttributes(desc, region, limit, "LFHCALLayerTiOVis");
  
  pvm = modScintAssembly.placeVolume(sliceti02H_vol, Position(0, 0, 0 ));
  pvm = modScintAssembly.placeVolume(sliceti02t_vol, Position(0, 0.5*(towerHeight+t_sep+h_nodge), 0 ));
  pvm = modScintAssembly.placeVolume(sliceti02b_vol, Position(0, -0.5*(towerHeight+t_sep+h_nodge), 0 ));

  // 8M module placement of scintillator for tower
  double rotZ[4] = {0,    0, 0,     0   };
  double rotY[4] = {M_PI, 0, M_PI,  0   };
  double rotX[4] = {0,    0, M_PI,  M_PI};
  double posX[4] = {(towerWidth*0.5+0.5*t_sep),     -(towerWidth*0.5+0.5*t_sep),  
                    (towerWidth*0.5+0.5*t_sep),     -(towerWidth*0.5+0.5*t_sep)};
  double posY[4] = {0.5*(towerHeight)+0.5*t_sep,    0.5*(towerHeight)+0.5*t_sep,
                    -(0.5*(towerHeight)+0.5*t_sep), -(0.5*(towerHeight)+0.5*t_sep)};
  double posZ[4] = {0,                              0,                  
                    0,                              0};
  // loop over all towers within same module
  for (int i = 0; i < 4; i++){
    Assembly modScintTowerAss = createScintillatorTower( desc,  basename+ _toString(i, "_tower_%d"),  
                                                            towerWidth, towerHeight, t_slice,
                                                            h_nodge, w_nodge, 
                                                            slice_mat, region, limit, vis, sens);
    pvm = modScintAssembly.placeVolume(modScintTowerAss, Transform3D(RotationZYX(rotZ[i], rotY[i], rotX[i]), Position(posX[i], posY[i], posZ[i] )));
    pvm.addPhysVolID("module", modID).addPhysVolID("tower", i).addPhysVolID("layer", layerID);
  }
  
  return modScintAssembly;         
}          


//************************************************************************************************************
//************************** create 8M module assembly  ******************************************************
//************************************************************************************************************
Assembly createEightMModule ( Detector& desc,
                              moduleParamsStrct mod_params,
                              std::vector<sliceParamsStrct> sl_params,
                              int modID,
                              double length,
                              SensitiveDetector sens
){
  std::string baseName = "LFHCAL_8M"+_toString(modID, "_%d");
  
  // assembly definition
  Assembly moduleAssembly(baseName);
  moduleAssembly.setVisAttributes(mod_params.mod_visStr);

  // placement operator
  PlacedVolume pvm;
  
  // ********************************************************************************
  // Casing definition
  // ********************************************************************************
  // geom definition 8M module casing
  Box         modFrontPlate( mod_params.mod_width / 2., mod_params.mod_height / 2., mod_params.mod_FWThick / 2.);
  Box         modSidePlateL( mod_params.mod_SWThick / 2., mod_params.mod_height / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);
  Box         modSidePlateR( mod_params.mod_SWThick / 2., mod_params.mod_height / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);
  Box         modTopPlate( (mod_params.mod_width-2*mod_params.mod_SWThick) / 2., mod_params.mod_TWThick / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);
  Box         modBottomPlate( (mod_params.mod_width-2*mod_params.mod_SWThick) / 2., mod_params.mod_TWThick / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);  
  Box         modBackCutOut( mod_params.mod_BIwidth / 2., mod_params.mod_BIheight / 2., mod_params.mod_BWThick / 2.);
  Box         modBackPlateFull( mod_params.mod_width / 2., mod_params.mod_height / 2., mod_params.mod_BWThick / 2.);
  SubtractionSolid modBackPlate(modBackPlateFull, modBackCutOut);
  
  // volume definition 8M module casing
  Volume  vol_modFrontPlate(baseName+"_FrontPlate",modFrontPlate,desc.material("Steel235"));
  Volume  vol_modBackPlate(baseName+"_BackPlate",modBackPlate,desc.material("Steel235"));
  Volume  vol_modSidePlateL(baseName+"_LeftSidePlate",modSidePlateL,desc.material("Steel235"));
  Volume  vol_modSidePlateR(baseName+"_RightSidePlate",modSidePlateR,desc.material("Steel235"));
  Volume  vol_modTopPlate(baseName+"_TopPlate",modTopPlate,desc.material("Steel235"));
  Volume  vol_modBottomPlate(baseName+"_BottomPlate",modBottomPlate,desc.material("Steel235"));
  
  vol_modFrontPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
  vol_modBackPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
  vol_modSidePlateL.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
  vol_modSidePlateR.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
  vol_modTopPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
  vol_modBottomPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
  
  int    layer_num = 1;
  
  double slice_z   = -length/2+mod_params.mod_FWThick; // Keeps track of layers' local z locations
  
  // Looping through the number of repeated layers & slices in each section
  for (int i = 0; i < (int)sl_params.size(); i++){
    slice_z += sl_params[i].slice_offset + sl_params[i].slice_thick / 2.; // Going to halfway point in layer

    //*************************************************
    // absorber plates
    //*************************************************
    Material slice_mat = desc.material(sl_params[i].slice_matStr);
    if (sl_params[i].slice_partID == 1 ){
      Assembly modAbsAssembly = createAbsorberPlateEightM( desc, 
                                                          baseName+"_AbsAssembly"+_toString(sl_params[i].layer_ID, "_layer_%d")+_toString(sl_params[i].slice_ID, "slice_%d"),
                                                          mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                          sl_params[i].slice_thick, mod_params.mod_nodgeDepth, mod_params.mod_nodgeWidthAbsA, mod_params.mod_nodgeWidthAbsB, mod_params.mod_nodgeWidthAbsC,
                                                          slice_mat, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr);
      // Placing slice within layer
      pvm = moduleAssembly.placeVolume(modAbsAssembly, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
    //*************************************************  
    // air & tyvek
    //*************************************************  
    } else if (sl_params[i].slice_partID == 2 ){
      Assembly modFillAssembly =  createFillerPlateEightM( desc, 
                                                          baseName+"_FillAssembly"+_toString(sl_params[i].layer_ID, "_layer_%d")+_toString(sl_params[i].slice_ID, "slice_%d"),
                                                          mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                          sl_params[i].slice_thick, mod_params.mod_nodgeDepth, mod_params.mod_nodgeWidthScin,
                                                          slice_mat, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr);
      // Placing slice within layer
      pvm = moduleAssembly.placeVolume(modFillAssembly, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
    //*************************************************
    // scintillator
    //*************************************************
    } else {
      Assembly modScintAssembly =  createScintillatorPlateEightM(  desc,
                                                                  baseName+"_ScintAssembly"+_toString(sl_params[i].layer_ID, "_layer_%d")+_toString(sl_params[i].slice_ID, "slice_%d"),
                                                                  modID,  layer_num, 
                                                                  mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                                  sl_params[i].slice_thick, mod_params.mod_nodgeDepth, mod_params.mod_nodgeWidthScin, mod_params.mod_sepDepth, 
                                                                  slice_mat, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr, sens);
      // Placing slice within layer
      pvm = moduleAssembly.placeVolume(modScintAssembly, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
    }
  }
  
  // placement 8M module casing
  pvm = moduleAssembly.placeVolume(vol_modFrontPlate, Position(0, 0, -( length-mod_params.mod_FWThick) / 2. ));
  pvm = moduleAssembly.placeVolume(vol_modBackPlate, Position(0, 0, ( length-mod_params.mod_BWThick) / 2. ));
  pvm = moduleAssembly.placeVolume(vol_modSidePlateL, Position(-(mod_params.mod_width-mod_params.mod_SWThick)/2., 0,  (mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  pvm = moduleAssembly.placeVolume(vol_modSidePlateR, Position((mod_params.mod_width-mod_params.mod_SWThick)/2., 0,(mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  pvm = moduleAssembly.placeVolume(vol_modTopPlate, Position(0, (mod_params.mod_height-mod_params.mod_TWThick)/2., (mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  pvm = moduleAssembly.placeVolume(vol_modBottomPlate, Position(0, -(mod_params.mod_height-mod_params.mod_TWThick)/2., (mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  
  return moduleAssembly;
}


//************************************************************************************************************
//************************** create 8M module assembly  ******************************************************
//************************************************************************************************************
Assembly createFourMModule ( Detector& desc,
                              moduleParamsStrct mod_params,
                              std::vector<sliceParamsStrct> sl_params,
                              int modID,
                              double length,
                              SensitiveDetector sens
){

  std::string baseName = "LFHCAL_4M"+_toString(modID, "_%d");
  
  // assembly definition
  Assembly moduleAssembly(baseName);
  moduleAssembly.setVisAttributes(mod_params.mod_visStr);

  // placement operator
  PlacedVolume pvm;
  
  // ********************************************************************************
  // Casing definition
  // ********************************************************************************
  // geom definition 8M module casing
  Box         modFrontPlate( mod_params.mod_width / 2., mod_params.mod_height / 2., mod_params.mod_FWThick / 2.);
  Box         modSidePlateL( mod_params.mod_SWThick / 2., mod_params.mod_height / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);
  Box         modSidePlateR( mod_params.mod_SWThick / 2., mod_params.mod_height / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);
  Box         modTopPlate( (mod_params.mod_width-2*mod_params.mod_SWThick) / 2., mod_params.mod_TWThick / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);
  Box         modBottomPlate( (mod_params.mod_width-2*mod_params.mod_SWThick) / 2., mod_params.mod_TWThick / 2., (length-mod_params.mod_FWThick-mod_params.mod_BWThick) / 2.);  
  Box         modBackCutOut( mod_params.mod_BIwidth / 2., mod_params.mod_BIheight / 2., mod_params.mod_BWThick / 2.);
  Box         modBackPlateFull( mod_params.mod_width / 2., mod_params.mod_height / 2., mod_params.mod_BWThick / 2.);
  SubtractionSolid modBackPlate(modBackPlateFull, modBackCutOut);
  
  // volume definition 8M module casing
  Volume  vol_modFrontPlate(baseName+"_FrontPlate",modFrontPlate,desc.material("Steel235"));
  Volume  vol_modBackPlate(baseName+"_BackPlate",modBackPlate,desc.material("Steel235"));
  Volume  vol_modSidePlateL(baseName+"_LeftSidePlate",modSidePlateL,desc.material("Steel235"));
  Volume  vol_modSidePlateR(baseName+"_RightSidePlate",modSidePlateR,desc.material("Steel235"));
  Volume  vol_modTopPlate(baseName+"_TopPlate",modTopPlate,desc.material("Steel235"));
  Volume  vol_modBottomPlate(baseName+"_BottomPlate",modBottomPlate,desc.material("Steel235"));
  
  vol_modFrontPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
  vol_modBackPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
  vol_modSidePlateL.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
  vol_modSidePlateR.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
  vol_modTopPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
  vol_modBottomPlate.setAttributes(desc, mod_params.mod_regStr, mod_params.mod_limStr, mod_params.mod_visStr);
  
  int    layer_num = 1;
  double slice_z   = -length/2+mod_params.mod_FWThick; // Keeps track of layers' local z locations
  
  // Looping through the number of repeated layers & slices in each section
  for (int i = 0; i < (int)sl_params.size(); i++){
    slice_z += sl_params[i].slice_offset + sl_params[i].slice_thick / 2.; // Going to halfway point in layer

    //*************************************************
    // absorber plates
    //*************************************************
    Material slice_mat = desc.material(sl_params[i].slice_matStr);
    if (sl_params[i].slice_partID == 1 ){
      Assembly modAbsAssembly = createAbsorberPlateFourM( desc, baseName+"_AbsAssembly"+_toString(sl_params[i].layer_ID, "_layer_%d")+_toString(sl_params[i].slice_ID, "slice_%d"),
                                                              mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                              sl_params[i].slice_thick, mod_params.mod_nodgeDepth, mod_params.mod_nodgeWidthAbsA, mod_params.mod_nodgeWidthAbsC,
                                                              slice_mat, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr);
      // Placing slice within layer
      pvm = moduleAssembly.placeVolume(modAbsAssembly, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
    //*************************************************  
    // air & tyvek
    //*************************************************  
    } else if (sl_params[i].slice_partID == 2 ){
      Assembly modFillAssembly =  createFillerPlateFourM( desc, baseName+"_FillAssembly"+_toString(sl_params[i].layer_ID, "_layer_%d")+_toString(sl_params[i].slice_ID, "slice_%d"),
                                                              mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                              sl_params[i].slice_thick, mod_params.mod_nodgeDepth, mod_params.mod_nodgeWidthScin,
                                                              slice_mat, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr);
      // Placing slice within layer
      pvm = moduleAssembly.placeVolume(modFillAssembly, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
    //*************************************************
    // scintillator
    //*************************************************
    } else {
      Assembly modScintAssembly =  createScintillatorPlateFourM( desc,baseName+"_ScintAssembly"+_toString(sl_params[i].layer_ID, "_layer_%d")+_toString(sl_params[i].slice_ID, "slice_%d"),
                                                                  modID,  layer_num, 
                                                                  mod_params.mod_height, mod_params.mod_width, mod_params.mod_TWThick, mod_params.mod_SWThick,
                                                                  sl_params[i].slice_thick, mod_params.mod_nodgeDepth, mod_params.mod_nodgeWidthScin, mod_params.mod_sepDepth, 
                                                                  slice_mat, sl_params[i].slice_regStr, sl_params[i].slice_limStr, sl_params[i].slice_visStr, sens);
      // Placing slice within layer
      pvm = moduleAssembly.placeVolume(modScintAssembly, Transform3D(RotationZYX(0, 0, 0), Position(0., 0., slice_z)));
    }
  }
  
  // placement 8M module casing
  pvm = moduleAssembly.placeVolume(vol_modFrontPlate, Position(0, 0, -( length-mod_params.mod_FWThick) / 2. ));
  pvm = moduleAssembly.placeVolume(vol_modBackPlate, Position(0, 0, ( length-mod_params.mod_BWThick) / 2. ));
  pvm = moduleAssembly.placeVolume(vol_modSidePlateL, Position(-(mod_params.mod_width-mod_params.mod_SWThick)/2., 0,  (mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  pvm = moduleAssembly.placeVolume(vol_modSidePlateR, Position((mod_params.mod_width-mod_params.mod_SWThick)/2., 0,(mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  pvm = moduleAssembly.placeVolume(vol_modTopPlate, Position(0, (mod_params.mod_height-mod_params.mod_TWThick)/2., (mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  pvm = moduleAssembly.placeVolume(vol_modBottomPlate, Position(0, -(mod_params.mod_height-mod_params.mod_TWThick)/2., (mod_params.mod_FWThick-mod_params.mod_BWThick)/2));
  
  return moduleAssembly;
}

//********************************************************************************************
//*                                                                                          *
//*                              Create detector                                             *
//==============================  MAIN FUNCTION  ============================================= 
//*                                                                                          *
//********************************************************************************************
static Ref_t createDetector(Detector& desc, xml_h handle, SensitiveDetector sens)
{
  // global detector variables
  xml_det_t   detElem = handle;
  std::string detName = detElem.nameStr();
  int         detID   = detElem.id();

  // general detector dimensions
  xml_dim_t dim    = detElem.dimensions();
  double    length = dim.z();    // Size along z-axis
  xml_dim_t pos = detElem.position();

  // 8M module specific loading
  xml_comp_t  eightM_xml        = detElem.child(_Unicode(eightmodule));
  xml_dim_t eightMmod_dim        = eightM_xml.dimensions();
  moduleParamsStrct eightM_params(getAttrOrDefault(eightMmod_dim, _Unicode(widthBackInner), 0.), 
                                  getAttrOrDefault(eightMmod_dim, _Unicode(heightBackInner), 0.), 
                                  getAttrOrDefault(eightMmod_dim, _Unicode(widthSideWall), 0.), 
                                  getAttrOrDefault(eightMmod_dim, _Unicode(widthTopWall), 0.), 
                                  getAttrOrDefault(eightMmod_dim, _Unicode(thicknessFrontWall), 0.), 
                                  getAttrOrDefault(eightMmod_dim, _Unicode(thicknessBackWall), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(width), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(height), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(nodgeWidthAbsA), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(nodgeWidthAbsB), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(nodgeWidthAbsC), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(nodgeWidthScin), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(nodgeDepth), 0.),
                                  getAttrOrDefault(eightMmod_dim, _Unicode(sepDepth), 0.), 
                                  eightM_xml.visStr(),
                                  eightM_xml.regionStr(),
                                  eightM_xml.limitsStr()
                                 );
  
  // 4M module specific loading
  xml_comp_t  fourM_xml        = detElem.child(_Unicode(fourmodule));
  xml_dim_t fourMmod_dim        = fourM_xml.dimensions();
  moduleParamsStrct fourM_params( getAttrOrDefault(fourMmod_dim, _Unicode(widthBackInner), 0.), 
                                  getAttrOrDefault(fourMmod_dim, _Unicode(heightBackInner), 0.), 
                                  getAttrOrDefault(fourMmod_dim, _Unicode(widthSideWall), 0.), 
                                  getAttrOrDefault(fourMmod_dim, _Unicode(widthTopWall), 0.), 
                                  getAttrOrDefault(fourMmod_dim, _Unicode(thicknessFrontWall), 0.), 
                                  getAttrOrDefault(fourMmod_dim, _Unicode(thicknessBackWall), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(width), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(height), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(nodgeWidthAbsA), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(nodgeWidthAbsB), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(nodgeWidthAbsC), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(nodgeWidthScin), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(nodgeDepth), 0.),
                                  getAttrOrDefault(fourMmod_dim, _Unicode(sepDepth), 0.), 
                                  fourM_xml.visStr(),
                                  fourM_xml.regionStr(),
                                  fourM_xml.limitsStr());

  std::vector<sliceParamsStrct> slice_Params;
  int layer_num = 1;
  for (xml_coll_t c(detElem, _U(layer)); c; ++c) {
    xml_comp_t x_layer         = c;
    int        repeat          = x_layer.repeat();
    // Looping through the number of repeated layers in each section
    for (int i = 0; i < repeat; i++) {
      int    slice_num = 1;
      
      // Looping over each layer's slices
      for (xml_coll_t l(x_layer, _U(slice)); l; ++l) {
        xml_comp_t  x_slice         = l;
        sliceParamsStrct slice_param( layer_num,  slice_num, getAttrOrDefault(l, _Unicode(type), 0.), x_slice.thickness(), getAttrOrDefault(l, _Unicode(offset), 0.), 
                                      x_slice.materialStr(), x_slice.visStr(), x_slice.regionStr(), x_slice.limitsStr());
        slice_Params.push_back(slice_param);
        ++slice_num;  
      }
    }
    layer_num++;
  }

  
  PlacedVolume phv;
  
  // create mother volume  
  DetElement det(detName, detID);
  Volume     motherVol = desc.pickMotherVolume(det);

  // create 8M modules
  int    moduleID   = 0;
  std::vector<double> xpos8M;
  std::vector<double> ypos8M;
  std::vector<double> zpos8M;
  
  for(xml_coll_t i(handle, _Unicode(eightmodulepositions)); i; ++i){
    xml_comp_t  x_mtrx = i;

    std::string   mtrx_name       = getAttrOrDefault<std::string>(x_mtrx, _Unicode(name), " ");
    std::string   mtrx_values     = getAttrOrDefault<std::string>(x_mtrx, _Unicode(values), " ");

    std::vector<double> *aptr = NULL;

    if(mtrx_name == "xpos")
      aptr = &xpos8M;
    else if(mtrx_name == "ypos")
      aptr = &ypos8M;
    else if(mtrx_name == "zpos")
      aptr = &zpos8M;
    else{
      printout(WARNING, "LFHCAL", "unknown <eightmodulepositions> data!");
      continue;
    }

    std::string delimiter = " ";
    size_t posC = 0;
    std::string token;
    while ((posC = mtrx_values.find(delimiter)) != std::string::npos) {
      token = mtrx_values.substr(0, posC);
      aptr->push_back(atof(token.c_str()));
      mtrx_values.erase(0, posC + delimiter.length());
    }
    aptr->push_back(atof(mtrx_values.c_str()));
  }

  if (xpos8M.size() != ypos8M.size() || xpos8M.size() != zpos8M.size()){
    std::cout << xpos8M.size() << "\t" <<ypos8M.size() <<  "\t" << zpos8M.size() << std::endl;
    std::cout <<  "idiot you can't count" << std::endl;
  } else {
    for (int e = 0; e < (int)xpos8M.size(); e++){
      std::cout <<  "LFHCAL placing 8M module: " << e << "/"<< (int)xpos8M.size() << "\t" << xpos8M[e] << "\t" << ypos8M[e] << "\t" << zpos8M[e]<< std::endl;
      Assembly  eightMassembly = createEightMModule ( desc, eightM_params, slice_Params, moduleID, length, sens);
      
      // Placing modules in world volume
      auto tr8M = Transform3D(Position(pos.x()-xpos8M[e]*dd4hep::mm-0.5*eightM_params.mod_width, pos.y()-ypos8M[e]*dd4hep::mm-0.5*eightM_params.mod_height, pos.z() +zpos8M[e]*dd4hep::mm + length / 2.));
      phv = motherVol.placeVolume(eightMassembly, tr8M);
      moduleID++;
    }
  }

  std::vector<double> xpos4M;
  std::vector<double> ypos4M;
  std::vector<double> zpos4M;

  for(xml_coll_t i(handle, _Unicode(fourmodulepositions)); i; ++i){
    xml_comp_t  x_mtrx = i;

    std::string   mtrx_name       = getAttrOrDefault<std::string>(x_mtrx, _Unicode(name), " ");
    std::string   mtrx_values     = getAttrOrDefault<std::string>(x_mtrx, _Unicode(values), " ");

    std::vector<double> *aptr = NULL;

    if(mtrx_name == "xpos")
      aptr = &xpos4M;
    else if(mtrx_name == "ypos")
      aptr = &ypos4M;
    else if(mtrx_name == "zpos")
      aptr = &zpos4M;
    else{
      printout(WARNING, "LFHCAL", "unknown <fourmodulepositions> data!");
      continue;
    }

    std::string delimiter = " ";
    size_t posC = 0;
    std::string token;
    while ((posC = mtrx_values.find(delimiter)) != std::string::npos) {
      token = mtrx_values.substr(0, posC);
      aptr->push_back(atof(token.c_str()));
      mtrx_values.erase(0, posC + delimiter.length());
    }
    aptr->push_back(atof(mtrx_values.c_str()));
  }
  
  // create 4M modules
  if (xpos4M.size() != ypos4M.size() || xpos4M.size() != zpos4M.size()){
    std::cout << xpos4M.size() << "\t" <<ypos4M.size() <<  "\t" << zpos4M.size() << std::endl;
    std::cout <<  "idiot you can't count" << std::endl;
  } else {
    for (int f = 0; f < (int)xpos4M.size(); f++){
      std::cout <<  "LFHCAL placing 4M module: " << f << "/"<< (int)xpos4M.size() << "\t" << xpos4M[f] << "\t" << ypos4M[f] << "\t" << zpos4M[f]<< std::endl;
      Assembly  fourMassembly = createFourMModule ( desc, fourM_params, slice_Params,  moduleID, length, sens);
    // Placing modules in world volume
      auto tr4M = Transform3D(Position(pos.x()-xpos4M[f]*dd4hep::mm-0.5*fourM_params.mod_width, pos.y()-ypos4M[f]*dd4hep::mm-0.5*fourM_params.mod_height, pos.z() +zpos4M[f]*dd4hep::mm + length / 2.));
      phv = motherVol.placeVolume(fourMassembly, tr4M);
      moduleID++;
    }
  }
  phv.addPhysVolID("system", detID);
  det.setPlacement(phv);

  return det;
}
DECLARE_DETELEMENT(epic_LFHCAL, createDetector)
