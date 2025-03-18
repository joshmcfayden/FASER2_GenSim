
setup_dict={

    "F2-default":{
        "name":"FASER2 Orig",# (Default)",
        #"color":"maroon",
        "color":"firebrick",
        "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
        "length":5,
        "distance":480,
        "channels": None,
        "G4":["FASER2_HepMC_v4_FASER2_Default_3rdTrkStation"],
        "effs":{"sep":[(0.1,'sep>0.1mm','abs(ep_y-em_y)>0.1'),
                       (1,'sep>1mm','abs(ep_y-em_y)>1'),
                       (5,'sep>5mm','abs(ep_y-em_y)>5'),
                       (10,'sep>10mm','abs(ep_y-em_y)>10'),
                      (100,'sep>100mm','abs(ep_y-em_y)>100')]}
    },

#    "S2-L10-D2":{
#        #"name":"S2 L=10m D=2m",
#        "name":"Old Baseline",
#        "color":"royalblue",
#        #"color":"darkgreen",
#        #"color":"orange",
#        "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
#        "length":10,
#        "distance":615,
#        "channels": None,
#        "G4":["FASER2_HepMC_v4_FASER2_Cavern_3rdTrkStation"],
#        "effs":{"sep":[(0.1,'sep>0.1mm','abs(ep_y-em_y)>0.1'),
#                       (1,'sep>1mm','abs(ep_y-em_y)>1'),
#                       (5,'sep>5mm','abs(ep_y-em_y)>5'),
#                       (10,'sep>10mm','abs(ep_y-em_y)>10'),
#                       (100,'sep>100mm','abs(ep_y-em_y)>100')]}
#
#    },
#
    
        "R1-L10-R1x3":{
            #"name":"R1 L=10m X=3m Y=1m",
            "name":"New Baseline (X=3m Y=1m)",
            #"color":"cornflowerblue",
            "color":"forestgreen",
            #"color":"limegreen",
            "style":"dashed",
            "selection":"(np.sqrt(x.x**2)<1.5) * (np.sqrt(x.y**2)<0.5)",
            "length":10,
            "distance":615,
            "channels": None,
             "G4":["FASER2_HepMC_v4_FASER2_Cavern_Rect_Baseline_Bhoriz_AllTrkStations"],
            "stations":['1','3'],
            "effs":{
#                "sep":[(0.1,'sep>0.1mm','abs(ep_x-em_x)>0.1'),
#                       (1,'sep>1mm','abs(ep_x-em_x)>1'),
#                       (5,'sep>5mm','abs(ep_x-em_x)>5'),
#                       (10,'sep>10mm','abs(ep_x-em_x)>10'),
#                       (100,'sep>100mm','abs(ep_x-em_x)>100')],
                "env":[(-250,'env<-250mm','abs(ep_x)<(1500-250)&&abs(em_x)<(1500-250)&&abs(ep_y)<(500-250)&&abs(em_y)<(500-250)'),
                       (-100,'env<-100mm','abs(ep_x)<(1500-100)&&abs(em_x)<(1500-100)&&abs(ep_y)<(500-100)&&abs(em_y)<(500-100)'),
                       (0,'env<0mm','abs(ep_x)<(1500)&&abs(em_x)<(1500)&&abs(ep_y)<(500)&&abs(em_y)<(500)'),
                       (100,'env<100mm','abs(ep_x)<(1500+100)&&abs(em_x)<(1500+100)&&abs(ep_y)<(500+100)&&abs(em_y)<(500+100)'),
                       (250,'env<250mm','abs(ep_x)<(1500+250)&&abs(em_x)<(1500+250)&&abs(ep_y)<(500+250)&&abs(em_y)<(500+250)'),
                       (500,'env<500mm','abs(ep_x)<(1500+500)&&abs(em_x)<(1500+500)&&abs(ep_y)<(500+500)&&abs(em_y)<(500+500)'),
                       (1000,'env<1000mm','abs(ep_x)<(1500+1000)&&abs(em_x)<(1500+1000)&&abs(ep_y)<(500+1000)&&abs(em_y)<(500+1000)')],
#                "envx":[(-1000,'env<-1000mm','abs(ep_x)<(1500-1000)&&abs(em_x)<(1500-1000)'),
#                        (-500,'env<-500mm','abs(ep_x)<(1500-500)&&abs(em_x)<(1500-500)'),
#                        (-250,'env<-250mm','abs(ep_x)<(1500-250)&&abs(em_x)<(1500-250)'),
#                        (-100,'env<-100mm','abs(ep_x)<(1500-100)&&abs(em_x)<(1500-100)'),
#                        (0,'env<0mm','abs(ep_x)<(1500)&&abs(em_x)<(1500)'),
#                        (100,'env<100mm','abs(ep_x)<(1500+100)&&abs(em_x)<(1500+100)'),
#                        (250,'env<250mm','abs(ep_x)<(1500+250)&&abs(em_x)<(1500+250)'),
#                        (500,'env<500mm','abs(ep_x)<(1500+500)&&abs(em_x)<(1500+500)'),
#                        (1000,'env<1000mm','abs(ep_x)<(1500+1000)&&abs(em_x)<(1500+1000)')],

            }
        },

#      "R1-L10-R0p5x3":{
#          #"name":"R1 L=10m X=3m Y=0.5m",
#          "name":"New Baseline (X=3m Y=0.5m)",
#          #"color":"lightsteelblue",
#          "color":"orchid",
#          #"color":"limegreen",
#          "style":"dotted",
#          "selection":"(np.sqrt(x.x**2)<1.5) * (np.sqrt(x.y**2)<0.25)",
#          "length":10,
#          "distance":615,
#          "channels": None,
#          "G4":["FASER2_HepMC_v4_FASER2_Cavern_Rect_KEKRect_3rdTrkStation"],
#          "effs":{"sep":[(0.1,'sep>0.1mm','abs(ep_x-em_x)>0.1'),
#                         (1,'sep>1mm','abs(ep_x-em_x)>1'),
#                         (5,'sep>5mm','abs(ep_x-em_x)>5'),
#                         (10,'sep>10mm','abs(ep_x-em_x)>10'),
#                         (100,'sep>100mm','abs(ep_x-em_x)>100')]}
#      },

#    "R1-L10-R0p5x2":{
#        #"name":"R1 L=10m X=2m Y=0.5m",
#        "name":"New Baseline (X=2m Y=0.5m)",
#        #"color":"lightsteelblue",
#        "color":"pink",
#        "style":"dotted",
#        "selection":"(np.sqrt(x.x**2)<1.0) * (np.sqrt(x.y**2)<0.25)",
#        "length":10,
#        "distance":615,
#        "channels": None,
#        "G4":["FASER2_HepMC_v4_FASER2_Cavern_Rect_KEKCircle_3rdTrkStation"],
#        "effs":{"sep":[(0.1,'sep>0.1mm','abs(ep_x-em_x)>0.1'),
#                       (1,'sep>1mm','abs(ep_x-em_x)>1'),
#                       (5,'sep>5mm','abs(ep_x-em_x)>5'),
#                       (10,'sep>10mm','abs(ep_x-em_x)>10'),
#                       (100,'sep>100mm','abs(ep_x-em_x)>100')]}
#    },
    

    
}
