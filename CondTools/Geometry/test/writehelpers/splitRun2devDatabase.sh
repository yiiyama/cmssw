#!/bin/sh
conddb_import -f sqlite_file:myfile.db -c sqlite_file:GeometryFileExtended2015dev.db -t XMLFILE_Geometry_TagXX_Extended2015dev_mc -i XMLFILE_Geometry_TagXX_Extended2015dev_mc
conddb_import -f sqlite_file:myfile.db -c sqlite_file:GeometryFileIdeal2015dev.db -t XMLFILE_Geometry_TagXX_Ideal2015dev_mc -i XMLFILE_Geometry_TagXX_Ideal2015dev_mc
conddb_import -f sqlite_file:myfile.db -c sqlite_file:TKRECO_Geometry.db -t TKRECO_Geometry2015dev_TagXX -i TKRECO_Geometry2015dev_TagXX
conddb_import -f sqlite_file:myfile.db -c sqlite_file:TKExtra_Geometry.db -t TKExtra_Geometry2015dev_TagXX -i TKExtra_Geometry2015dev_TagXX
conddb_import -f sqlite_file:myfile.db -c sqlite_file:TKParameters_Geometry.db -t TKParameters_Geometry2015dev_TagXX -i TKParameters_Geometry2015dev_TagXX
conddb_import -f sqlite_file:myfile.db -c sqlite_file:EBRECO_Geometry.db -t EBRECO_Geometry2015dev_TagXX -i EBRECO_Geometry2015dev_TagXX
conddb_import -f sqlite_file:myfile.db -c sqlite_file:EERECO_Geometry.db -t EERECO_Geometry2015dev_TagXX -i EERECO_Geometry2015dev_TagXX
conddb_import -f sqlite_file:myfile.db -c sqlite_file:EPRECO_Geometry.db -t EPRECO_Geometry2015dev_TagXX -i EPRECO_Geometry2015dev_TagXX
conddb_import -f sqlite_file:myfile.db -c sqlite_file:HCALRECO_Geometry.db -t HCALRECO_Geometry2015dev_TagXX -i HCALRECO_Geometry2015dev_TagXX
conddb_import -f sqlite_file:myfile.db -c sqlite_file:HCALParameters_Geometry.db -t HCALParameters_Geometry2015dev_TagXX -i HCALParameters_Geometry2015dev_TagXX
conddb_import -f sqlite_file:myfile.db -c sqlite_file:CTRECO_Geometry.db -t CTRECO_Geometry2015dev_TagXX -i CTRECO_Geometry2015dev_TagXX
conddb_import -f sqlite_file:myfile.db -c sqlite_file:ZDCRECO_Geometry.db -t ZDCRECO_Geometry2015dev_TagXX -i ZDCRECO_Geometry2015dev_TagXX
conddb_import -f sqlite_file:myfile.db -c sqlite_file:CASTORRECO_Geometry.db -t CASTORRECO_Geometry2015dev_TagXX -i CASTORRECO_Geometry2015dev_TagXX
conddb_import -f sqlite_file:myfile.db -c sqlite_file:CSCRECO_Geometry.db -t CSCRECO_Geometry2015dev_TagXX -i CSCRECO_Geometry2015dev_TagXX
conddb_import -f sqlite_file:myfile.db -c sqlite_file:CSCRECODIGI_Geometry.db -t CSCRECODIGI_Geometry2015dev_TagXX -i CSCRECODIGI_Geometry2015dev_TagXX
conddb_import -f sqlite_file:myfile.db -c sqlite_file:DTRECO_Geometry.db -t DTRECO_Geometry2015dev_TagXX -i DTRECO_Geometry2015dev_TagXX
conddb_import -f sqlite_file:myfile.db -c sqlite_file:RPCRECO_Geometry.db -t RPCRECO_Geometry2015dev_TagXX -i RPCRECO_Geometry2015dev_TagXX
