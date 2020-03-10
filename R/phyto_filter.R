


phyto_filter <- function(flowfile, pig_channels, com_channels, ph = 0.05, width_channel = "SSC.W") {

  #gate based on the pigments
  pigment_gating <- pigment_gate(flowfile = flowfile, pig_channels = pig_channels,
                                 ph = ph, width_channel = width_channel)

  #gate based on the complexity channels
  sizecomp_gating <- complexity_gate(flowfile = flowfile, com_channels = com_channels,
                                     ph = ph, width_channel = width_channel
                                    )

  ### formulate a master full flow file

 }
