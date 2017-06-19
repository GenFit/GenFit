#include <ConstField.h>
#include <TGeoMaterialInterface.h>
#include <FieldManager.h>
#include <MaterialEffects.h>

#include <TGeoManager.h>


int main() {

    // init geometry and mag. field
    new TGeoManager("Geometry", "Geane geometry");
    TGeoManager::Import("genfitGeom.root");
    genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0., 1.5));
    genfit::FieldManager::getInstance()->useCache(true, 8);
    genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

    genfit::MaterialEffects::getInstance()->drawdEdx(11);

    return 0;
}