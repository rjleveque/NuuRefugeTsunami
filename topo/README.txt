
------------------------------------------------
Run
    python fetch_topo_nuu19.py
to fetch 1/9 arcsecond tile and crop it to the desired region.
Note that the tile 
    ncei19_n19x00_w155x75_2021v1.tif
is fetched from:
    https://chs.coast.noaa.gov/htdata/raster2/elevation/NCEI_ninth_Topobathy_Hawaii_9428/
and is 159 MB.


------------------------------------------------
Run
    python fetch_hawaii6sec.py
to fetch a subset of the Hawaii 6 second DEM described at
    https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.ngdc.mgg.dem:460/html
    
------------------------------------------------
Run
    python make_topo_plots.py
to create a plot showing the topo and contours of the 10m slump.
