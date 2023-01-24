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
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include <XML/Helper.h>
using namespace dd4hep;

//************************************************************************************************************
//************************** create 8M module assembly  ******************************************************
//************************************************************************************************************
Volume createEightMModule(Detector& desc,
                          int modID, double length, SensitiveDetector sens)
{
  double modBox_length_tot = length;
  double modBox_length = modBox_length_tot-10*cm;
  double modBox_steel_length = 120*cm;
  double modBox_tungsten_length = 10*cm;
  double modBox_width = 20*cm;
  double modBox_height = 10*cm;
  
  std::string baseName = "GFHCAL_8M" + _toString(modID, "_mod_%d");
  PlacedVolume pvm;

  int visDetails = 0;
  visDetails = 1;
  // visDetails = 2;

  // std::cout << "create 8M module" << std::endl;
  // std::cout << "length = " << length << std::endl;
  Box    modBox(modBox_width / 2., modBox_height / 2., length / 2.);
  Volume vol_modBox(baseName, modBox, desc.material("Air"));
  if(visDetails){
    vol_modBox.setVisAttributes(desc.visAttributes("InvisibleWithDaughters"));
  }else{
    vol_modBox.setVisAttributes(desc.visAttributes("AnlViolet"));
  }

  double tyvek_thickness = 0.34*mm;
  double sciSeg_length_tot = 5*cm;
  double sciSeg_width_tot = 4*mm+tyvek_thickness*2;
  double sciSeg_height_tot = 5*cm;

  double sciSeg_length = 5*cm-tyvek_thickness*2;
  double sciSeg_width = 4*mm;
  double sciSeg_height = 5*cm-tyvek_thickness*2;


  double modBox_sidewall_thickness = 2 * mm;
  // double modBox_topwall_thickness = 0.5* mm;
  double absorber_thickness = 14.0* mm;
  double pcb_thickness = 1.0* mm;
  double miniframe_thickness = 1.0*mm;


  Box    scintBox(sciSeg_width / 2., sciSeg_height / 2., sciSeg_length / 2.);
  Volume vol_scintBox(baseName + "_Scintillator8M" , scintBox, desc.material("Polystyrene"));
  if(visDetails==2){
    vol_scintBox.setVisAttributes(desc.visAttributes("GFHCALLayerScintVis"));
  }else{
    vol_scintBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }
  // vol_scintBox.setSensitiveDetector(sens);

  Box    steelWallBox(modBox_sidewall_thickness / 2., modBox_height / 2., modBox_length_tot / 2.);
  Volume vol_steelWallBox(baseName + "_FeWall8M" , steelWallBox, desc.material("Steel235"));
  if(visDetails){
    vol_steelWallBox.setVisAttributes(desc.visAttributes("AnlOrange"));
  }else{
    vol_steelWallBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    steelAbsorberBox(absorber_thickness / 2., modBox_height / 2., modBox_steel_length / 2.);
  Volume vol_steelAbsorberBox(baseName + "_FeAbsorber8M" , steelAbsorberBox, desc.material("Steel235"));
  if(visDetails){
    vol_steelAbsorberBox.setVisAttributes(desc.visAttributes("AnlLight_Gray"));
  }else{
    vol_steelAbsorberBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    tungstenAbsorberBox(absorber_thickness / 2., modBox_height / 2., modBox_tungsten_length / 2.);
  Volume vol_tungstenAbsorberBox(baseName + "_WAbsorber8M" , tungstenAbsorberBox, desc.material("Tungsten"));
  if(visDetails){
    vol_tungstenAbsorberBox.setVisAttributes(desc.visAttributes("AnlViolet"));
  }else{
    vol_tungstenAbsorberBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    steelMiniFrameBox(miniframe_thickness / 2., modBox_height / 2., modBox_length_tot / 2.);
  Volume vol_steelMiniFrameBox(baseName + "_MiniFrame8M" , steelMiniFrameBox, desc.material("Steel235"));
  if(visDetails==2){
    vol_steelMiniFrameBox.setVisAttributes(desc.visAttributes("AnlDarkRed"));
  }else{
    vol_steelMiniFrameBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    pcbBox(pcb_thickness / 2., modBox_height / 2., modBox_length_tot / 2.);
  Volume vol_pcbBox(baseName + "_PCB8M" , pcbBox, desc.material("Fr4"));
  if(visDetails==2){
    vol_pcbBox.setVisAttributes(desc.visAttributes("AnlDarkGreen"));
  }else{
    vol_pcbBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    tyvekBox(tyvek_thickness / 2., modBox_height / 2., modBox_length / 2.);
  Volume vol_tyvekBox(baseName + "_Tyvek8M" , tyvekBox, desc.material("Tyvek"));
  if(visDetails==2){
    vol_tyvekBox.setVisAttributes(desc.visAttributes("AnlGold"));
  }else{
    vol_tyvekBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  double miniBox_thickness = absorber_thickness+pcb_thickness+sciSeg_width_tot+miniframe_thickness*2;
  int nLayers_x = (int)( (modBox_width-2*modBox_sidewall_thickness) / miniBox_thickness );
  int   nLayers_z      = (int)(modBox_length/sciSeg_length_tot);
  int   nLayers_y      = 2;

  vol_modBox.placeVolume(vol_steelWallBox, Transform3D(RotationZYX(0, 0, 0),
    Position(-modBox_width/2 + modBox_sidewall_thickness/2, 0, 0)));
  vol_modBox.placeVolume(vol_steelWallBox, Transform3D(RotationZYX(0, 0, 0),
    Position(modBox_width/2 - modBox_sidewall_thickness/2, 0, 0)));

  for(int ilx = 0; ilx < nLayers_x; ilx++)
  {
    vol_modBox.placeVolume(vol_steelAbsorberBox, Transform3D(RotationZYX(0, 0, 0),
      Position(-modBox_width/2 + modBox_sidewall_thickness + absorber_thickness/2 + ilx*miniBox_thickness, 0, 0)));

    vol_modBox.placeVolume(vol_tungstenAbsorberBox, Transform3D(RotationZYX(0, 0, 0),
      Position(-modBox_width/2 + modBox_sidewall_thickness + absorber_thickness/2 + ilx*miniBox_thickness, 0, -modBox_length_tot/2+modBox_tungsten_length/2)));

    vol_modBox.placeVolume(vol_pcbBox, Transform3D(RotationZYX(0, 0, 0),
      Position(-modBox_width/2 + modBox_sidewall_thickness + absorber_thickness + miniframe_thickness + pcb_thickness/2 + ilx*miniBox_thickness, 0, 0)));

    vol_modBox.placeVolume(vol_steelMiniFrameBox, Transform3D(RotationZYX(0, 0, 0),
      Position(-modBox_width/2 + modBox_sidewall_thickness + absorber_thickness + miniframe_thickness/2 + ilx*miniBox_thickness, 0, 0)));
    vol_modBox.placeVolume(vol_steelMiniFrameBox, Transform3D(RotationZYX(0, 0, 0),
      Position(-modBox_width/2 + modBox_sidewall_thickness + miniBox_thickness - miniframe_thickness/2 + ilx*miniBox_thickness, 0, 0)));

    vol_modBox.placeVolume(vol_tyvekBox, Transform3D(RotationZYX(0, 0, 0),
      Position(-modBox_width/2 + modBox_sidewall_thickness + absorber_thickness + miniframe_thickness + pcb_thickness + tyvek_thickness/2 + ilx*miniBox_thickness, 0, -modBox_tungsten_length/2)));
    vol_modBox.placeVolume(vol_tyvekBox, Transform3D(RotationZYX(0, 0, 0),
      Position(-modBox_width/2 + modBox_sidewall_thickness + miniBox_thickness - miniframe_thickness - tyvek_thickness/2 + ilx*miniBox_thickness, 0, -modBox_tungsten_length/2)));
    for(int ily = 0; ily < nLayers_y; ily++)
    {
      for(int ilz = 0; ilz < nLayers_z; ilz++)
      {
        // Volume vol_scintBox(baseName + _toString(ilx, "_ix_%d") + _toString(ily, "_iy_%d") + _toString(ilz, "_iz_%d") , scintBox, desc.material("Polystyrene"));
        // if(visDetails){
        //   vol_scintBox.setVisAttributes(desc.visAttributes("GFHCALLayerScintVis"));
        // }else{
        //   vol_scintBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
        // }
        // vol_scintBox.setSensitiveDetector(sens);

        sens.setType("calorimeter");
        vol_scintBox.setSensitiveDetector(sens);
        pvm = vol_modBox.placeVolume(vol_scintBox, Transform3D(RotationZYX(0, 0, 0), 
        Position(
          -modBox_width/2 + modBox_sidewall_thickness + absorber_thickness + miniframe_thickness + pcb_thickness + sciSeg_width_tot/2  + ilx*miniBox_thickness, 
          -modBox_height/2 + sciSeg_height_tot/2 + ily*sciSeg_height_tot,
          -modBox_length_tot/2 + sciSeg_length_tot/2 + ilz*sciSeg_length_tot)));
        pvm.addPhysVolID("module", modID);
        pvm.addPhysVolID("layerx", ilx);
        pvm.addPhysVolID("layery", ily);
        pvm.addPhysVolID("layerz", ilz);
      }
    }
  }


  return vol_modBox;
}

//************************************************************************************************************
//************************** create 4M module assembly  ******************************************************
//************************************************************************************************************
Volume createFourMModule(Detector& desc,
                          int modID, double length, SensitiveDetector sens)
{
  double modBox_length_tot = length;
  double modBox_length = modBox_length_tot-10*cm;
  double modBox_steel_length = 120*cm;
  double modBox_tungsten_length = 10*cm;
  double modBox_width = 10*cm;
  double modBox_height = 10*cm;
  
  std::string baseName = "GFHCAL_4M" + _toString(modID, "_mod_%d");
  PlacedVolume pvm;


  int visDetails = 0;
  visDetails = 1;
  // int visDetails = 2;

  // std::cout << "create 4M module" << std::endl;
  // std::cout << "length = " << length << std::endl;
  Box    modBox(modBox_width / 2., modBox_height / 2., length / 2.);
  Volume vol_modBox(baseName, modBox, desc.material("Air"));
  if(visDetails){
    vol_modBox.setVisAttributes(desc.visAttributes("InvisibleWithDaughters"));
  } else {
    vol_modBox.setVisAttributes(desc.visAttributes("AnlOrange"));
  }
  double tyvek_thickness = 0.34*mm;
  double sciSeg_length_tot = 5*cm;
  double sciSeg_width_tot = 4*mm+tyvek_thickness*2;
  double sciSeg_height_tot = 5*cm;

  double sciSeg_length = 5*cm-tyvek_thickness*2;
  double sciSeg_width = 4*mm;
  double sciSeg_height = 5*cm-tyvek_thickness*2;


  double modBox_sidewall_thickness = 2 * mm;
  // double modBox_topwall_thickness = 0.5* mm;
  double absorber_thickness = 14.0* mm;
  double pcb_thickness = 1.0* mm;
  double miniframe_thickness = 1.0*mm;


  Box    scintBox(sciSeg_width / 2., sciSeg_height / 2., sciSeg_length / 2.);
  // Volume vol_scintBox(baseName + _toString(ilx, "_ix_%d") + _toString(ily, "_iy_%d") + _toString(ilz, "_iz_%d") , scintBox, desc.material("Polystyrene"));
  Volume vol_scintBox(baseName + "_Scintillator4M" , scintBox, desc.material("Polystyrene"));
  if(visDetails==2){
    vol_scintBox.setVisAttributes(desc.visAttributes("GFHCALLayerScintVis"));
  }else{
    vol_scintBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    steelWallBox(modBox_sidewall_thickness / 2., modBox_height / 2., modBox_length_tot / 2.);
  Volume vol_steelWallBox(baseName + "_FeWall4M" , steelWallBox, desc.material("Steel235"));
  if(visDetails){
    vol_steelWallBox.setVisAttributes(desc.visAttributes("AnlOrange"));
  }else{
    vol_steelWallBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    steelAbsorberBox(absorber_thickness / 2., modBox_height / 2., modBox_steel_length / 2.);
  Volume vol_steelAbsorberBox(baseName + "_FeAbsorber4M" , steelAbsorberBox, desc.material("Steel235"));
  if(visDetails){
    vol_steelAbsorberBox.setVisAttributes(desc.visAttributes("AnlLight_Gray"));
  }else{
    vol_steelAbsorberBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    tungstenAbsorberBox(absorber_thickness / 2., modBox_height / 2., modBox_tungsten_length / 2.);
  Volume vol_tungstenAbsorberBox(baseName + "_WAbsorber4M" , tungstenAbsorberBox, desc.material("Tungsten"));
  if(visDetails){
    vol_tungstenAbsorberBox.setVisAttributes(desc.visAttributes("AnlBlue"));
  }else{
    vol_tungstenAbsorberBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    steelMiniFrameBox(miniframe_thickness / 2., modBox_height / 2., modBox_length_tot / 2.);
  Volume vol_steelMiniFrameBox(baseName + "_MiniFrame4M" , steelMiniFrameBox, desc.material("Steel235"));
  if(visDetails==2){
    vol_steelMiniFrameBox.setVisAttributes(desc.visAttributes("AnlDarkRed"));
  }else{
    vol_steelMiniFrameBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    pcbBox(pcb_thickness / 2., modBox_height / 2., modBox_length_tot / 2.);
  Volume vol_pcbBox(baseName + "_PCB4M" , pcbBox, desc.material("Fr4"));
  if(visDetails==2){
    vol_pcbBox.setVisAttributes(desc.visAttributes("AnlDarkGreen"));
  }else{
    vol_pcbBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  Box    tyvekBox(tyvek_thickness / 2., modBox_height / 2., modBox_length / 2.);
  Volume vol_tyvekBox(baseName + "_Tyvek4M" , tyvekBox, desc.material("Tyvek"));
  if(visDetails==2){
    vol_tyvekBox.setVisAttributes(desc.visAttributes("AnlGold"));
  }else{
    vol_tyvekBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  double miniBox_thickness = absorber_thickness+pcb_thickness+sciSeg_width_tot+miniframe_thickness*2;
  int nLayers_x = (int)( (modBox_width-2*modBox_sidewall_thickness) / miniBox_thickness );
  int   nLayers_z      = (int)(modBox_length/sciSeg_length_tot);
  int   nLayers_y      = 2;

  double addAbsorber_thickness = (modBox_width-2*modBox_sidewall_thickness) - nLayers_x*miniBox_thickness;
  Box    steelAbsorberBoxAdd(addAbsorber_thickness / 2., modBox_height / 2., modBox_steel_length / 2.);
  Volume vol_steelAbsorberBoxAdd(baseName + "_FeAbsorberAdd4M" , steelAbsorberBoxAdd, desc.material("Steel235"));
  if(visDetails){
    vol_steelAbsorberBoxAdd.setVisAttributes(desc.visAttributes("AnlBrown"));
  }else{
    vol_steelAbsorberBoxAdd.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }
  Box    tungstenAbsorberBoxAdd(addAbsorber_thickness / 2., modBox_height / 2., modBox_tungsten_length / 2.);
  Volume vol_tungstenAbsorberBoxAdd(baseName + "_WAbsorberAdd4M" , tungstenAbsorberBoxAdd, desc.material("Tungsten"));
  if(visDetails){
    vol_tungstenAbsorberBoxAdd.setVisAttributes(desc.visAttributes("AnlTeal"));
  }else{
    vol_tungstenAbsorberBoxAdd.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
  }

  vol_modBox.placeVolume(vol_steelWallBox, Transform3D(RotationZYX(0, 0, 0),
    Position(-modBox_width/2 + modBox_sidewall_thickness/2, 0, 0)));
  vol_modBox.placeVolume(vol_steelWallBox, Transform3D(RotationZYX(0, 0, 0),
    Position(modBox_width/2 - modBox_sidewall_thickness/2, 0, 0)));

  vol_modBox.placeVolume(vol_steelAbsorberBoxAdd, Transform3D(RotationZYX(0, 0, 0),
    Position(modBox_width/2 - modBox_sidewall_thickness -addAbsorber_thickness/2, 0, 0)));
  vol_modBox.placeVolume(vol_tungstenAbsorberBoxAdd, Transform3D(RotationZYX(0, 0, 0),
    Position(modBox_width/2 - modBox_sidewall_thickness -addAbsorber_thickness/2, 0, -modBox_length_tot/2+modBox_tungsten_length/2)));
  

  for(int ilx = 0; ilx < nLayers_x; ilx++)
  {
    vol_modBox.placeVolume(vol_steelAbsorberBox, Transform3D(RotationZYX(0, 0, 0),
      Position(-modBox_width/2 + modBox_sidewall_thickness + absorber_thickness/2 + ilx*miniBox_thickness, 0, 0)));

    vol_modBox.placeVolume(vol_tungstenAbsorberBox, Transform3D(RotationZYX(0, 0, 0),
      Position(-modBox_width/2 + modBox_sidewall_thickness + absorber_thickness/2 + ilx*miniBox_thickness, 0, -modBox_length_tot/2+modBox_tungsten_length/2)));

    vol_modBox.placeVolume(vol_pcbBox, Transform3D(RotationZYX(0, 0, 0),
      Position(-modBox_width/2 + modBox_sidewall_thickness + absorber_thickness + miniframe_thickness + pcb_thickness/2 + ilx*miniBox_thickness, 0, 0)));

    vol_modBox.placeVolume(vol_steelMiniFrameBox, Transform3D(RotationZYX(0, 0, 0),
      Position(-modBox_width/2 + modBox_sidewall_thickness + absorber_thickness + miniframe_thickness/2 + ilx*miniBox_thickness, 0, 0)));
    vol_modBox.placeVolume(vol_steelMiniFrameBox, Transform3D(RotationZYX(0, 0, 0),
      Position(-modBox_width/2 + modBox_sidewall_thickness + miniBox_thickness - miniframe_thickness/2 + ilx*miniBox_thickness, 0, 0)));

    vol_modBox.placeVolume(vol_tyvekBox, Transform3D(RotationZYX(0, 0, 0),
      Position(-modBox_width/2 + modBox_sidewall_thickness + absorber_thickness + miniframe_thickness + pcb_thickness + tyvek_thickness/2 + ilx*miniBox_thickness, 0, -modBox_tungsten_length/2)));
    vol_modBox.placeVolume(vol_tyvekBox, Transform3D(RotationZYX(0, 0, 0),
      Position(-modBox_width/2 + modBox_sidewall_thickness + miniBox_thickness - miniframe_thickness - tyvek_thickness/2 + ilx*miniBox_thickness, 0, -modBox_tungsten_length/2)));
    for(int ily = 0; ily < nLayers_y; ily++)
    {
      for(int ilz = 0; ilz < nLayers_z; ilz++)
      {
        // Volume vol_scintBox(baseName + _toString(ilx, "_ix_%d") + _toString(ily, "_iy_%d") + _toString(ilz, "_iz_%d") , scintBox, desc.material("Polystyrene"));
        // if(visDetails){
        //   vol_scintBox.setVisAttributes(desc.visAttributes("GFHCALLayerScintVis"));
        // }else{
        //   vol_scintBox.setVisAttributes(desc.visAttributes("InvisibleNoDaughters"));
        // }
        // vol_scintBox.setSensitiveDetector(sens);

        sens.setType("calorimeter");
        vol_scintBox.setSensitiveDetector(sens);

        pvm = vol_modBox.placeVolume(vol_scintBox, Transform3D(RotationZYX(0, 0, 0), 
        Position(
          -modBox_width/2 + modBox_sidewall_thickness + absorber_thickness + miniframe_thickness + pcb_thickness + sciSeg_width_tot/2  + ilx*miniBox_thickness, 
          -modBox_height/2 + sciSeg_height_tot/2 + ily*sciSeg_height_tot,
          -modBox_length_tot/2 + sciSeg_length_tot/2 + ilz*sciSeg_length_tot)));
        pvm.addPhysVolID("module", modID);
        pvm.addPhysVolID("layerx", ilx);
        pvm.addPhysVolID("layery", ily);
        pvm.addPhysVolID("layerz", ilz);
      }
    }
  }


  return vol_modBox;
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
  xml_det_t   x_det = handle;
  int       det_id   = x_det.id();
  std::string    det_name = x_det.nameStr();

  xml_dim_t dim    = x_det.dimensions();
  double    length = dim.z(); // Size along z-axis
  xml_dim_t pos    = x_det.position();

  std::cout << "global GFHCAL position" << pos.x() << "\t" << pos.y() << "\t" << pos.z() << std::endl;




  DetElement sdet(det_name, det_id);

  Assembly assembly(det_name);
  PlacedVolume phv;

  int moduleID = 0;



  // phv = motherVol.placeVolume(eightMassembly, Position(0, 0, 0));


  // create 8M modules
  std::vector<double> xpos8M;
  std::vector<double> ypos8M;
  std::vector<double> zpos8M;

  for (xml_coll_t i(handle, _Unicode(eightmodulepositions)); i; ++i) {
    xml_comp_t x_mtrx = i;

    std::string mtrx_name   = getAttrOrDefault<std::string>(x_mtrx, _Unicode(name), " ");
    std::string mtrx_values = getAttrOrDefault<std::string>(x_mtrx, _Unicode(values), " ");

    std::vector<double>* aptr = NULL;

    if (mtrx_name == "xpos")
      aptr = &xpos8M;
    else if (mtrx_name == "ypos")
      aptr = &ypos8M;
    else if (mtrx_name == "zpos")
      aptr = &zpos8M;
    else {
      printout(WARNING, "GFHCAL", "unknown <eightmodulepositions> data!");
      continue;
    }

    std::string delimiter = " ";
    size_t      posC      = 0;
    std::string token;
    while ((posC = mtrx_values.find(delimiter)) != std::string::npos) {
      token = mtrx_values.substr(0, posC);
      aptr->push_back(atof(token.c_str()));
      mtrx_values.erase(0, posC + delimiter.length());
    }
    aptr->push_back(atof(mtrx_values.c_str()));
  }

  if (xpos8M.size() != ypos8M.size() || xpos8M.size() != zpos8M.size()) {
    std::cout << xpos8M.size() << "\t" << ypos8M.size() << "\t" << zpos8M.size() << std::endl;
    std::cout << "idiot you can't count" << std::endl;
  } else {
    for (int e = 0; e < (int)xpos8M.size(); e++) {
      if (e % 20 == 0)
        std::cout << "GFHCAL placing 8M module: " << e << "/" << (int)xpos8M.size() << "\t" << xpos8M[e] << "\t"
                  << ypos8M[e] << "\t" << zpos8M[e] << std::endl;

      Volume eightMassembly =
                createEightMModule(desc, moduleID, length,sens);




      auto tr8M =
          Transform3D(Position(pos.x() - xpos8M[e] * dd4hep::cm - 0.5 * 20*cm,
                               pos.y() - ypos8M[e] * dd4hep::cm, pos.z() + zpos8M[e] * dd4hep::cm + length / 2.));
      phv = assembly.placeVolume(eightMassembly, tr8M);
      phv.addPhysVolID("module", moduleID);

      moduleID++;
    }
  }

  // create 4M modules
  std::vector<double> xpos4M;
  std::vector<double> ypos4M;
  std::vector<double> zpos4M;

  for (xml_coll_t i(handle, _Unicode(fourmodulepositions)); i; ++i) {
    xml_comp_t x_mtrx = i;

    std::string mtrx_name   = getAttrOrDefault<std::string>(x_mtrx, _Unicode(name), " ");
    std::string mtrx_values = getAttrOrDefault<std::string>(x_mtrx, _Unicode(values), " ");

    std::vector<double>* aptr = NULL;

    if (mtrx_name == "xpos")
      aptr = &xpos4M;
    else if (mtrx_name == "ypos")
      aptr = &ypos4M;
    else if (mtrx_name == "zpos")
      aptr = &zpos4M;
    else {
      printout(WARNING, "GFHCAL", "unknown <fourmodulepositions> data!");
      continue;
    }

    std::string delimiter = " ";
    size_t      posC      = 0;
    std::string token;
    while ((posC = mtrx_values.find(delimiter)) != std::string::npos) {
      token = mtrx_values.substr(0, posC);
      aptr->push_back(atof(token.c_str()));
      mtrx_values.erase(0, posC + delimiter.length());
    }
    aptr->push_back(atof(mtrx_values.c_str()));
  }

  if (xpos4M.size() != ypos4M.size() || xpos4M.size() != zpos4M.size()) {
    std::cout << xpos4M.size() << "\t" << ypos4M.size() << "\t" << zpos4M.size() << std::endl;
    std::cout << "idiot you can't count" << std::endl;
  } else {
    for (int e = 0; e < (int)xpos4M.size(); e++) {
      if (e % 20 == 0)
        std::cout << "GFHCAL placing 4M module: " << e << "/" << (int)xpos4M.size() << "\t" << xpos4M[e] << "\t"
                  << ypos4M[e] << "\t" << zpos4M[e] << std::endl;

      Volume eightMassembly =
                createFourMModule(desc, moduleID, length,sens);

      auto tr4M =
          Transform3D(Position(pos.x() - xpos4M[e] * dd4hep::cm - 0.5 * 10*cm,
                               pos.y() - ypos4M[e] * dd4hep::cm, pos.z() + zpos4M[e] * dd4hep::cm + length / 2.));
      phv = assembly.placeVolume(eightMassembly, tr4M);
      phv.addPhysVolID("module", moduleID);

      moduleID++;
    }
  }

  Volume     motherVol = desc.pickMotherVolume(sdet);
  // std::cout << motherVol.name() << std::endl;

  Transform3D  tr        = Translation3D(0., 0., 0.) * RotationZYX(0.,0.,0.);
  PlacedVolume envPV     = motherVol.placeVolume(assembly, tr);
  envPV.addPhysVolID("system", det_id);
  sdet.setPlacement(envPV);

  return sdet;
}
DECLARE_DETELEMENT(epic_GFHCAL, createDetector)
