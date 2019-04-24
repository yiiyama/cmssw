#include "L1Trigger/L1THGCal/interface/veryfrontend/HGCalVFESummationImpl.h"

HGCalVFESummationImpl::
HGCalVFESummationImpl(const edm::ParameterSet& conf):
  thickness_corrections_(conf.getParameter<std::vector<double>>("ThicknessCorrections")),
  LSB_silicon_fC_(conf.getParameter<double>("siliconCellLSB_fC")),
  LSB_scintillator_MIP_(conf.getParameter<double>("scintillatorCellLSB_MIP")),
  thresholds_silicon_(conf.getParameter<std::vector<double>>("thresholdsSilicon")),
  threshold_scintillator_(conf.getParameter<double>("thresholdScintillator"))
{
  if(thickness_corrections_.size()!=3)
  {
    throw cms::Exception("Configuration") <<
      thickness_corrections_.size() << " thickness corrections are given instead of 3 (the number of sensor thicknesses)";
  }
  if(thresholds_silicon_.size()!=3)
  {
    throw cms::Exception("Configuration") <<
      thresholds_silicon_.size() << " silicon thresholds are given instead of 3 (the number of sensor thicknesses)";
  }
}

void 
HGCalVFESummationImpl::
triggerCellSums(const HGCalTriggerGeometryBase& geometry, 
                const std::vector<std::pair<DetId, uint32_t > >& linearized_dataframes,
                std::unordered_map<uint32_t, uint32_t>& payload)
{
  if(linearized_dataframes.empty()) return;
  // sum energies in trigger cells
  for(const auto& frame : linearized_dataframes)
  {
    DetId cellid(frame.first);
    uint32_t value = frame.second;

    // Apply noise threshold before summing into trigger cells
    if(triggerTools_.isSilicon(cellid))
    {
      int thickness = triggerTools_.thicknessIndex(cellid);
      double threshold = thresholds_silicon_.at(thickness);
      value = ( value*LSB_silicon_fC_ > threshold ? value : 0 );
    }
    else if(triggerTools_.isScintillator(cellid))
    {
      value = ( value*LSB_scintillator_MIP_ > threshold_scintillator_ ? value : 0 );
    }
    if(value==0) continue;

    // find trigger cell associated to cell
    uint32_t tcid = geometry.getTriggerCellFromCell(cellid);
    payload.emplace(tcid, 0); // do nothing if key exists already

    // equalize value among cell thicknesses for Silicon parts
    if(triggerTools_.isSilicon(cellid))
    {
      int thickness = triggerTools_.thicknessIndex(cellid);
      double thickness_correction = thickness_corrections_.at(thickness);
      value = (double)value*thickness_correction;
    }

    // sums energy for the same trigger cell id
    payload[tcid] += value; // 32 bits integer should be largely enough 
  }

}
