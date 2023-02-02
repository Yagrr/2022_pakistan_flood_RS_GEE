/*
  ===========================================================================================
                       2022 PAKISTAN FLOOD MAPPING AND AFFECTED CROPLAND ESTIMATION
  ===========================================================================================
  
The following script is used for mapping the flood extent of the 2022 Pakistan flood
during the August month. The Sentinel-1 GRD C-band data product was used with a change
detection approach to detect flooding. In order to detect the affected cropland, the flood
extent map is combined with the Copernicus Global Land Cover map which is up to date till 2019.
  
Material:
  - Sentinel-1 GRD C-band. Preprocessed: Thermal-Noise removal, Radiometric calibration, Terrain-correction, speckle filter
        Timeframe:
      - Before flood - 1st March to 31st April
      - During flood - 1st August to 31st August
        
  - Copernicus Global Land Cover Layers: CGLS-LC100 Collection 3
        Timeframe: August 2019 
        
  - Administrative boundary of Pakistan shapefiles
        Share links
          Pakistan boundary:
          https://code.earthengine.google.com/?asset=projects/sirtas-gee/assets/0_pakistan
          
          Provincial boundary - Sindh Province, Area of Interest (AOI):
          https://code.earthengine.google.com/?asset=projects/sirtas-gee/assets/2_1_sindh_district_Karachi_merged
          
          Project Files Github repository link (.shp available):
          https://github.com/Yagrr/2022_pakistan_flood_RS_GEE
          
  NOTE: The Edge Otsu algorithm takes 4-5 minutes to run for the first time
  
  Index:
  1. Data Import
  2. Sentinel-1 Pre-Processing
    2.2 Sentinel-1 Post-Processing
    2.3 Flood statistics
  3. Export
  4.  Display Products
  5. Functions used in this script
  ===========================================================================================
                                          1. Data Import
  =========================================================================================== 
*/
// User input date and polarization

var beforeStart = '2022-05-01';
var beforeEnd = '2022-06-30';
var afterStart = '2022-08-01';
var afterEnd = '2022-08-31';

var polarization = 'VV';

// CGLC Accuracy threhsold
// Classification Probability must be >50%
    var accuracy_threshold = 50; 

// Import boundaries

var pakistan = pakistan_border, 
    sindh = sindh_border, // AOI
    pakistan_geometry = pakistan.geometry(),
    sindh_geometry = sindh.geometry();

Map.addLayer(pakistan,{color:'gray'},'Pakistan',false);  
Map.addLayer(sindh,{color:'gray'},'Sindh',false);
Map.centerObject(sindh,8);

/*
  ===========================================================================================
                                    2. Sentinel-1 Pre-Processing
  ===========================================================================================                                 
*/

// Load Sentinel-1 C-band S1 Ground Range collection (log scaling, VH cross-polar)

      var S1 = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(ee.Filter.eq('instrumentMode','IW'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', polarization))
        .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
        .filter(ee.Filter.eq('resolution_meters',10))
        .filterBounds(sindh) // Clip to AOI
        .select(polarization);
        
// Filter by date
      //Before flood: 1 March - 30 April 
      //During flood: 1-31 August
      var S1_dry_raw = S1.filterDate(beforeStart, beforeEnd);
      var S1_flood_raw = S1.filterDate(afterStart, afterEnd);
      
// S-1 Print metadata (Number of Scenes)
      // Print dates of before images to console
        var S1_dry_count = S1_dry_raw.size();
        print(ee.String('Number of Before Flood scenes: ').cat('(').cat(S1_dry_count).cat(')'),
          getdate(S1_dry_raw), S1_dry_raw);
        
      // Print dates of after images to console
        var after_count = S1_flood_raw.size();
        print(ee.String('Number of After Flood scenes: ').cat('(').cat(after_count).cat(')'),
          getdate(S1_flood_raw), S1_flood_raw);
          
      // Function to extract date from metadata
        function getdate(imagecol){
          var time_period = imagecol.reduceColumns(ee.Reducer.minMax(), ["system:time_start"]);
          var printed = ee.String('from ')
            .cat(ee.Date(time_period.get('min')).format('YYYY-MM-dd'))
            .cat(' to ')
            .cat(ee.Date(time_period.get('max')).format('YYYY-MM-dd'));
          return printed;
        }
      
// S-1 mosaic, clip by AOI and preview
     var S1_dry = S1_dry_raw.mosaic().clip(sindh),  
     S1_flood = S1_flood_raw.mosaic().clip(sindh);

    // Map.addLayer(S1_dry,{min:-25,max:-15}, 'Before Floods',false);
    // Map.addLayer(S1_flood,{min:-25,max:-15}, 'After Floods',false);
      
// Refined Lee Filter
    var S1_dry_filtered = ee.Image(powerToDb(refinedLee(dbToPower(S1_dry))));
    var S1_flood_filtered = ee.Image(powerToDb(refinedLee(dbToPower(S1_flood))));
    
// SAR composite data for Before and During flood images
        var S1_collectionVH = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(ee.Filter.eq('instrumentMode','IW'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VH'))
        .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING')) 
        .filter(ee.Filter.eq('resolution_meters',10))
        .filterBounds(sindh) // Clip to AOI
        .select('VH');
        
        var S1_collectionVV = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(ee.Filter.eq('instrumentMode','IW'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
        .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING')) 
        .filter(ee.Filter.eq('resolution_meters',10))
        .filterBounds(sindh) // Clip to AOI
        .select('VV');
      
      var S1_dry_VV = S1_collectionVV.filterDate(beforeStart, beforeEnd).mosaic().clip(sindh),
          S1_dry_VH = S1_collectionVH.filterDate(beforeStart, beforeEnd).mosaic().clip(sindh);
          
      var S1_flood_VV = S1_collectionVV.filterDate(afterStart, afterEnd).mosaic().clip(sindh),
          S1_flood_VH = S1_collectionVH.filterDate(afterStart, afterEnd).mosaic().clip(sindh);

// Filter
      var VV_dry = ee.Image(powerToDb(refinedLee(dbToPower(S1_dry_VV)))),
          VH_dry = ee.Image(powerToDb(refinedLee(dbToPower(S1_dry_VH))));
      
      var VV_flood = ee.Image(powerToDb(refinedLee(dbToPower(S1_flood_VV)))),
          VH_flood = ee.Image(powerToDb(refinedLee(dbToPower(S1_flood_VH))));
        
Map.addLayer(VV_dry.addBands(VH_dry).addBands(VV_dry), {min:-13,max:-7},'Before-flood SAR composite',false);
Map.addLayer(VV_flood.addBands(VH_flood).addBands(VV_flood), {min:-13,max:-7},'During flood SAR composite',false);
//Map.addLayer(VV_flood.addBands(VH_flood).addBands(VH_flood.divide(VV_flood) ), {min:-13,max:-7},'During flood SAR composite');
/*
  ===========================================================================================
                                   2.1 Sentinel-1 Flood Mapping
  ===========================================================================================                                 
*/

// Difference layer of Before and During flood for change detection
    var S1_diff = S1_flood_filtered.subtract(S1_dry_filtered);

    // Edge Otsu: Input keyword arguments (kwarg)
    var kwargDefaults = {
        'initialThreshold': 1, // Initial threshold set to 1
        'reductionScale': 180,
        'smoothing': 100,
        'bandName': 'constant',
        'connectedPixels': 30,
        'edgeLength': 20,
        'smoothEdges': 20,
        'cannyThreshold': 1,
        'cannySigma': 1,
        'cannyLt': 0.05,
        'maxBuckets': 255,
        'minBucketWidth': 0.001,
        'maxRaw': 1e6,
        'invert': false,
        'verbose': true
    };
    
// Acquiring Canny Edge detection, finding optimal threshold and mapping
    var flood_raw = edgeOtsu(S1_diff,kwargDefaults,sindh_geometry); 

/*
  ===========================================================================================
                                2.2 Sentinel-1 Post-Processing
  ===========================================================================================                                 
*/

// Eliminate false-positives

  // JRC Global Surface Water Mapping Layers
  // Permanent water mask
    var water = ee.Image('JRC/GSW1_0/GlobalSurfaceWater').select('seasonality');
    var datawater = water.gte(8).updateMask(water.gte(8)); // Water > 8 months
    
      // Binary Mask
      //Flooded layer where perennial water bodies are assigned a 0 value
      var water_mask = flood_raw.where(datawater,0);
      // final flooded area without pixels in perennial waterbodies
      var flood = water_mask.updateMask(water_mask);

  // Eliminate Noise - Filter pixels with 8 neighbours
        var connections = flood.connectedPixelCount();    
        var flood = flood.updateMask(connections.gte(8));
      
  // Mask out areas with more than 5% slope using a Digital Elevation Model 
      var DEM = ee.Image('WWF/HydroSHEDS/03VFDEM');
      var terrain = ee.Algorithms.Terrain(DEM);
      var slope = terrain.select('slope');
      
// Final flood extent
      var flood = flood.updateMask(slope.lt(5));
    

//  Cropland Mask

  // Copernicus Global Land Cover
  //  Discrete Classification layer and clip to Sindh province
    var CGLC = ee.Image("COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019") // Discrete Classification dataset
                 .select('discrete_classification').clip(sindh);
    var CGLC_proba = ee.Image("COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019") // Classification Quality
           .select('discrete_classification-proba').clip(sindh);
             
    // Accuracy threshold
      var CGLC_proba90 = CGLC_proba.gt(accuracy_threshold); 
      
  // Cropland mask. Must be classified as cropland with > 90 probability
      var crop_mask_raw = CGLC.eq(40).updateMask(CGLC_proba90); // Cropland = 40 (See discrete_classification Table in docs)
      var cropland = crop_mask_raw.updateMask(crop_mask_raw);  // Filter out 0 values
      
      var flooded_cropland = flood.updateMask(cropland);
      


// Precipitation Map
  
  var precipitation_raw = ee.ImageCollection("NASA/GPM_L3/IMERG_V06")
                          .filterBounds(sindh)
                          .filterDate(afterStart, afterEnd)
                          .select('precipitationCal');
                          
    // Add rainfall accumulation into map
    var rainfall_total= precipitation_raw.reduce(ee.Reducer.sum());
    var rainfall = rainfall_total.clip(sindh);

    // Rainfall vis parameter
    var precip_palette = [
      '000096','0064ff', '00b4ff', '33db80', '9beb4a',
      'ffeb00', 'ffb300', 'ff6400', 'eb1e00', 'af0000'
    ];
    
    var precipVis = {min: 0.0, max: 2000};
    Map.addLayer(rainfall, precipVis, "Sum of rainfalls over August", false);
    

// Final Product
    // SAR 
    Map.addLayer(S1_dry_filtered,{min:-15,max:-1}, 'Before Floods | Filtered',false); // Before flood
    Map.addLayer(S1_flood_filtered,{min:-15,max:-1}, 'After Floods | Filtered',false);// During flood
    Map.addLayer(S1_diff,{min:-8,max:0}, 'Difference | Filtered',false); // Difference
    // Classification
    Map.addLayer(cropland,{palette:"ForestGreen"},"Cropland");
    Map.addLayer(flood,{palette:"RoyalBlue"},"Flood Extent");
    Map.addLayer(flooded_cropland,{palette:"Crimson"},"Inundated cropland");

/*
  ===========================================================================================
                                      2.3 Flood statistics
  ===========================================================================================                                 
*/

// Get pixel area of flood extent
    var flood_pixelarea = flood
      .multiply(ee.Image.pixelArea()); //calculate the area of each pixel (m^2)
    
    var flood_stats = flood_pixelarea.reduceRegion({
      reducer: ee.Reducer.sum(), //sum all pixels with area information                
      geometry: sindh,
      scale: 250,
      maxPixels: 1e14
      });
    
    // convert area to hectares
    var flood_ha = flood_stats
      .getNumber("constant")
      .divide(10000)
      .round();

// Get pixel area of cropland
    var cropland_pixelarea = cropland
      .multiply(ee.Image.pixelArea());
    
    var cropland_stats = cropland_pixelarea.reduceRegion({
      reducer: ee.Reducer.sum(), //sum all pixels with area information                
      geometry: sindh,
      scale: 250,
      maxPixels: 1e14
      });
    
    // convert area to hectares
    var cropland_ha = cropland_stats
      .getNumber("discrete_classification")
      .divide(10000)
      .round();


// Get pixel area of affected cropland
    var floodcrop_pixelarea = flooded_cropland
      .multiply(ee.Image.pixelArea()); 
    
    // sum pixels of affected cropland layer
    var floodcrop_stats = floodcrop_pixelarea.reduceRegion({
      reducer: ee.Reducer.sum(), //sum all pixels with area information                
      geometry: sindh,
      scale: 250,
      maxPixels: 1e14
      });
    
    // convert area to hectares
    var floodcrop_ha = floodcrop_stats
      .getNumber("constant")
      .divide(10000)
      .round();
 
// Affected crop (%)
    var floodcrop_percent = floodcrop_ha.divide(cropland_ha).multiply(100).round();

// Convert to integer then string for easier manipulation
var str_flood_ha = flood_ha.toInt().format("%s"),
    str_cropland_ha = cropland_ha.toInt().format("%s"),
    str_floodcrop_ha = floodcrop_ha.toInt().format("%s"),
    str_floodcrop_percent = floodcrop_percent.toInt().format("%s");

// Print to console
print("Flood extent (ha): ",str_flood_ha); 
print("Cropland (ha): ",str_cropland_ha); 
print("Flooded cropland (ha): ",str_floodcrop_ha);
print("Percentage (%) of cropland flooded: ",str_floodcrop_percent);


/* 
  ===========================================================================================
                                          3. Export
  ===========================================================================================                                 
*/

// Before Flood
  Export.image.toDrive({image:S1_dry_filtered.visualize({min:-15,max:-1}), 
    description: 'before_flood',
    fileNamePrefix: 'before_flood',
    region: sindh, 
    scale:250,
    folder:'GEE_export',
    crs: 'EPSG:4326', 
    maxPixels: 1e13
  });
// During Flood
  Export.image.toDrive({image:S1_flood_filtered.visualize({min:-15,max:-1}), 
    description: 'during_flood',
    fileNamePrefix: 'during_flood',
    region: sindh, 
    scale:250,
    folder:'GEE_export',
    crs: 'EPSG:4326', 
    maxPixels: 1e13
  });
// Difference layer
  Export.image.toDrive({image:S1_diff.visualize({min:-8,max:0}), 
    description: 'difference',
    fileNamePrefix: 'difference',
    region: sindh, 
    scale:250,
    folder:'GEE_export',
    crs: 'EPSG:4326', 
    maxPixels: 1e13
  });

// // FCC Composite Before Flood
//   Export.image.toDrive({image:S1_dry_filtered.visualize({min:-15,max:-1}), 
//     description: 'before_flood',
//     fileNamePrefix: 'before_flood',
//     region: sindh, 
//     scale:250,
//     folder:'GEE_export',
//     crs: 'EPSG:4326', 
//     maxPixels: 1e13
//   });
// // FCC Composite During Flood
//   Export.image.toDrive({image:S1_dry_filtered.visualize({min:-15,max:-1}), 
//     description: 'before_flood',
//     fileNamePrefix: 'before_flood',
//     region: sindh, 
//     scale:250,
//     folder:'GEE_export',
//     crs: 'EPSG:4326', 
//     maxPixels: 1e13
//   });

// Flood extent
  Export.image.toDrive({
    image: flood.visualize({palette:"RoyalBlue"}), 
    description: 'Flood_extent',
    fileNamePrefix: 'flood_extent',
    region: sindh, 
    scale:250,
    folder:'GEE_export',
    crs: 'EPSG:4326', 
    maxPixels: 1e13
  });

// Cropland
  Export.image.toDrive({
    image: cropland.visualize({palette:"ForestGreen"}), 
    description: 'Agricultural_land',
    fileNamePrefix: 'cropland',
    region: sindh, 
    scale:250,
    folder:'GEE_export',
    crs: 'EPSG:4326', 
    maxPixels: 1e13
  });

// Flooded Cropland
  Export.image.toDrive({image: flooded_cropland.visualize({palette:"Crimson"}), 
    description: 'Flooded_agricultural_land',
    fileNamePrefix: 'flooded_cropland',
    region: sindh, 
    scale:250,
    folder:'GEE_export',
    crs: 'EPSG:4326', 
    maxPixels: 1e13
  });
  
//  Precipitation Map
  Export.image.toDrive({image: rainfall.select('precipitationCal_sum'), 
    description: 'rainfall_sum',
    fileNamePrefix: 'august_rainfall_sum',
    region: sindh, 
    scale:1000,
    folder:'GEE_export',
    crs: 'EPSG:4326', 
    maxPixels: 1e13
  });


/* 
  ===========================================================================================
                                     4.  Display Products
  ===========================================================================================                                 
*/
// Code derived from UN-Spider: https://code.earthengine.google.com/f5c2f984c053c8ea574bfcd4040d084e
// https://un-spider.org/advisory-support/recommended-practices/recommended-practice-google-earth-engine-flood-mapping/step-by-step

// Set position of panel where the results will be displayed 
var results = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px',
    width: '250px'
  }
});

//Prepare the visualization parameters of the labels 
var textVis = {
  'margin':'0px 8px 2px 0px',
  'fontWeight':'bold',
  };
var numberVIS = {
  'margin':'0px 0px 15px 0px', 
  'color':'Crimson',
  'fontWeight':'bold'
  };
var subTextVis = {
  'margin':'0px 0px 2px 0px',
  'fontSize':'12px',
  'color':'grey'
  };

var titleTextVis = {
  'margin':'0px 0px 15px 0px',
  'fontSize': '18px', 
  'font-weight':'bold', 
  'color': 'black'
  };

// Create lables of the results 
// Title and time period
var title = ui.Label('Results', titleTextVis);
var text1 = ui.Label('Flood status between:',textVis);
var number1 = ui.Label(afterStart.concat(" and ",afterEnd),numberVIS);

// Estimated flood extent 
var text2 = ui.Label('Estimated flood extent:',textVis);
var text2_1 =  ui.Label('based on Sentinel-1 Imagery',subTextVis);
var number2 = ui.Label('Please wait...',numberVIS); 
str_flood_ha.evaluate(function(val){number2.setValue(val+' hectares')}),numberVIS;

// Estimated area of affected cropland 
var text3 = ui.Label('Estimated affected cropland:',textVis);
var text3_1 =  ui.Label('based on Copernicus Global Land Cover (100m)',subTextVis);
var number3 = ui.Label('Please wait...',numberVIS);
str_floodcrop_ha.evaluate(function(val){number3.setValue(val+' hectares')}),numberVIS;
var number4 = ui.Label('Please wait...',numberVIS);
str_floodcrop_percent.evaluate(function(val){number4.setValue(val+'% of cropland affected')}),numberVIS;

// Data sourced from:
var textcredit = ui.Label('UI script derived from: UN-SPIDER December 2019', subTextVis);

// Add the labels to the panel 
results.add(ui.Panel([
        title,
        text1,
        number1,
        text2,
        text2_1,
        number2,
        text3,
        text3_1,
        number3,
        number4,
        textcredit
        ]
      ));

// Add the panel to the map 
Map.add(results);

// Displaying Legend

// Create legend (*credits to thisearthsite on Open Geo Blog: https://mygeoblog.com/2016/12/09/add-a-legend-to-to-your-gee-map/)
// set position of panel
var legend = ui.Panel({
  style: {
    position: 'bottom-right',
    padding: '8px 15px',
  }
});
 
// Create legend title
var legendTitle = ui.Label('Legend',titleTextVis);
 
// Add the title to the panel
legend.add(legendTitle);
 
// Creates and styles 1 row of the legend.
var makeRow = function(color, name) {
 
      // Create the label that is actually the colored box.
      var colorBox = ui.Label({
        style: {
          backgroundColor: color,
          // Use padding to give the box height and width.
          padding: '8px',
          margin: '0 0 4px 0'
        }
      });
 
      // Create the label filled with the description text.
      var description = ui.Label({
        value: name,
        style: {margin: '0 0 4px 6px'}
      });
 
      // return the panel
      return ui.Panel({
        widgets: [colorBox, description],
        layout: ui.Panel.Layout.Flow('horizontal')
      });
};
 
//  Palette with the colors
var palette =['RoyalBlue', 'ForestGreen', 'Crimson'];
 
// name of the legend
var names = ['Potentially flooded areas','Cropland','Affected Cropland'];
 
// Add color and and names
for (var i = 0; i < 3; i++) {
  legend.add(makeRow(palette[i], names[i]));
  }  

Map.add(legend);

/* 
  ===========================================================================================
                                5. Functions used in this script
  ===========================================================================================                                 
*/

/*
  ==========================
  Refined Lee Speckle Filter
  ==========================                           
*/

// Applying a Refined Lee Speckle filter as coded in the SNAP 3.0 S1TBX:
// https://github.com/senbox-org/s1tbx/blob/master/s1tbx-op-sar-processing/src/main/java/org/esa/s1tbx/sar/gpf/filtering/SpeckleFilters/RefinedLee.java
// Adapted to GEE by Guido Lemoine 
// Link to original code: https://code.earthengine.google.com/2ef38463ebaf5ae133a478f173fd0ab5

function powerToDb(img){
  return ee.Image(10).multiply(img.log10());
}

function dbToPower(img){
  return ee.Image(10).pow(img.divide(10));
}

// Refined Lee filter (3x3 kernel) adapted to GEE Code Editor by Guido Lemoine
 function refinedLee(image) {
  
  var bandNames = image.bandNames();
  image = dbToPower(image);
  
  var result = ee.ImageCollection(bandNames.map(function(b){
    var img = image.select([b]);
    
    // img must be in natural units, i.e. not in dB!
    // Set up 3x3 kernels 
    var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
    var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);
  
    var mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
    var variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);
  
    // Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
    var sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);
  
    var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);
  
    // Calculate mean and variance for the sampled windows and store as 9 bands
    var sample_mean = mean3.neighborhoodToBands(sample_kernel); 
    var sample_var = variance3.neighborhoodToBands(sample_kernel);
  
    // Determine the 4 gradients for the sampled windows
    var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
    gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
    gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
    gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());
  
    // And find the maximum gradient amongst gradient bands
    var max_gradient = gradients.reduce(ee.Reducer.max());
  
    // Create a mask for band pixels that are the maximum gradient
    var gradmask = gradients.eq(max_gradient);
  
    // duplicate gradmask bands: each gradient represents 2 directions
    gradmask = gradmask.addBands(gradmask);
  
    // Determine the 8 directions
    var directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1);
    directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2));
    directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3));
    directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4));
    // The next 4 are the not() of the previous 4
    directions = directions.addBands(directions.select(0).not().multiply(5));
    directions = directions.addBands(directions.select(1).not().multiply(6));
    directions = directions.addBands(directions.select(2).not().multiply(7));
    directions = directions.addBands(directions.select(3).not().multiply(8));
  
    // Mask all values that are not 1-8
    directions = directions.updateMask(gradmask);
  
    // "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
    directions = directions.reduce(ee.Reducer.sum());  
  
    //var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000'];
    //Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false);
  
    var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));
  
    // Calculate localNoiseVariance
    var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);
  
    // Set up the 7*7 kernels for directional statistics
    var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));
  
    var diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], 
      [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);
  
    var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
    var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);
  
    // Create stacks for mean and variance using the original kernels. Mask with relevant direction.
    var dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
    var dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));
  
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));
  
    // and add the bands for rotated kernels
    for (var i=1; i<4; i++) {
      dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
      dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
      dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
      dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    }
  
    // "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
    dir_mean = dir_mean.reduce(ee.Reducer.sum());
    dir_var = dir_var.reduce(ee.Reducer.sum());
  
    // A finally generate the filtered value
    var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));
  
    var K = varX.divide(dir_var);
  
    return dir_mean.add(K.multiply(img.subtract(dir_mean)))
      .arrayProject([0])
      // Get a multi-band image bands.
      .arrayFlatten([['sum']])
      .float();
  })).toBands().rename(bandNames);
  return powerToDb(ee.Image(result));
}


/*
  ===========================================
  Automatic thresholding with Edge Otsu method
  ===========================================                          
*/

// Edge Otsu adapted from Markert et al. (2020) DOI: https://doi.org/10.3390/rs12152469
// https://mygeoblog.com/2021/01/25/edge-otsu-for-surface-water-mapping-detection/
// Edge Otsu algorithm: Canny Edge + Otsu. (Dependency: Histogram Constructor, Otsu algo)

function edgeOtsu(imageCollection,kwargs,boundary) {
  
  // force imageCollection type for input, used for some subsequent processing
  imageCollection = ee.ImageCollection(imageCollection);
  // var geom = imageCollection.map(function(img){
  //   return img.geometry();
  // }).union(1).geometry();
  // print(geom);
  var geom = boundary;
  
  // get list of band names used later
  var bandList = ee.Image(imageCollection.first()).bandNames();
  
  
    // Default Edge Otsu inputs
  var kwargDefaults = {
        'initialThreshold': 1, // Initial threshold set to 1
        'reductionScale': 180,
        'smoothing': 100,
        'bandName': 'constant',
        'connectedPixels': 30,
        'edgeLength': 20,
        'smoothEdges': 20,
        'cannyThreshold': 1,
        'cannySigma': 1,
        'cannyLt': 0.05,
        'maxBuckets': 255,
        'minBucketWidth': 0.001,
        'maxRaw': 1e6,
        'invert': false,
        'verbose': true
    };

  // define default parameterization for keywords
  var kwargKeys = [];
  for(var key in kwargDefaults) kwargKeys.push( key );
  var params;
  var i,k,v;
  // loop through the keywords and construct ee.Dictionary from them,
  // if the key is defined in the input then pass else use default
  params = ee.Dictionary(kwargs);
  for (i=0;i<kwargKeys.length;i++) {
    k = kwargKeys[i];
    v = kwargDefaults[k];
    params = ee.Dictionary(
      ee.Algorithms.If(params.contains(k),params,params.set(k,v))
    );
  }
  
  // parameters for all methods
  var initialThreshold = ee.Number( params.get('initialThreshold') ),
      reductionScale   = ee.Number( params.get('reductionScale') ),
      smoothing        = ee.Number( params.get('smoothing') ),
      bandName         = ee.String( params.get('bandName') ),
      connectedPixels  = ee.Number( params.get('connectedPixels') ),
      edgeLength       = ee.Number( params.get('edgeLength') ),
      smoothEdges      = ee.Number( params.get('smoothEdges') ),
      cannyThreshold   = ee.Number( params.get('cannyThreshold') ),
      cannySigma       = ee.Number( params.get('cannySigma') ),
      cannyLt          = ee.Number( params.get('cannyLt') ),
      maxBuckets       = ee.Number( params.get('maxBuckets') ),
      minBucketWidth   = ee.Number( params.get('minBucketWidth') ),
      maxRaw           = ee.Number( params.get('maxRaw') ),
      invert           = params.get('invert'),
      verbose          = params.get('verbose').getInfo();
      

  var img = imageCollection.select(bandName).mean();
  
  // get preliminary water
  var binary = img.lt(initialThreshold).rename('binary');

  // get canny edges
  var canny = ee.Algorithms.CannyEdgeDetector(binary, cannyThreshold, cannySigma);
  // process canny edges
  var connected  = canny.updateMask(canny).lt(cannyLt).connectedPixelCount(connectedPixels, true);
  var edges      = connected.gte(edgeLength);
  edges          = edges.updateMask(edges);
  var edgeBuffer = edges.focal_max(smoothEdges, 'square', 'meters');

  // get histogram for Otsu
  var histogram_image = img.updateMask(edgeBuffer);

  // histogram_image = histogram_image.clip(geometry2)

  var histogram = ee.Dictionary(histogram_image.reduceRegion({
    reducer:ee.Reducer.histogram(maxBuckets, minBucketWidth,maxRaw)
      .combine('mean', null, true).combine('variance', null,true),
    geometry: geom,
    scale: reductionScale,
    maxPixels: 1e13,
    tileScale:16
  }).get(bandName.cat('_histogram')));
  
// Otsu algorithm
function otsu(histogram) {
  var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
  var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);
  
  var indices = ee.List.sequence(1, size);
  
  // Compute between sum of squares, where each mean partitions the data.
  var bss = indices.map(function(i) {
    var aCounts = counts.slice(0, 0, i);
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMeans = means.slice(0, 0, i);
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
    var bCount = total.subtract(aCount);
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
           bCount.multiply(bMean.subtract(mean).pow(2)));
  });
  // print(ui.Chart.array.values(ee.Array(bss), 0, means));
  // Return the mean value corresponding to the maximum BSS.
  return means.sort(bss).get([-1]);
}
  
  var threshold = otsu(histogram);
  
  // Histogram Constructor
function constructHistChart(histogram,threshold){
  var counts = ee.List(histogram.get('histogram'));
  var buckets = ee.List(histogram.get('bucketMeans'));
  // construct array for visualization of threshold in chart
  var segment = ee.List.repeat(0, counts.size());
  var maxFrequency   = ee.Number(counts.reduce(ee.Reducer.max()));
  var threshIndex    = buckets.indexOf(threshold);
  segment            = segment.set(threshIndex, maxFrequency);
  var histChart = ui.Chart.array.values(ee.Array.cat([counts, segment], 1), 0, buckets)
    .setSeriesNames(['Values', 'Threshold'])
    .setChartType('ColumnChart');
  return histChart;
}
  
  if (verbose){
    var chart = constructHistChart(histogram,threshold)
      .setOptions({
        title: 'Edge Search Histogram',
        hAxis: {
          title: 'Values',
        },
        vAxis:{
          title:'Count'
        } 
      });
    print('Algorithm parameters:',params);
    print("Calculated threshold:",threshold);
    print('Thresholding histogram:',chart);
  }
  
  // segment image and mask 0 values (not water)
  var waterImg = ee.Image(ee.Algorithms.If(invert,img.gt(threshold),img.lt(threshold)));
  //waterImg = waterImg.updateMask(waterImg); // Binary mask
  return waterImg;
}


/*
  ===========================================================================================
                                          (Unused)
                                      Sentinel-2 Visuals
  ===========================================================================================                                 
*/

//     var S2_SRdry= ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
//       .filterDate(beforeStart, beforeEnd)  //During flood: March-April 
//       .filterBounds(pakistan);

//     var S2_SRflood = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
//       .filterDate(afterStart, afterEnd)  //During flood: 1-31 August
//       .filterBounds(pakistan);
  
// // // Cloud Mask: https://developers.google.com/earth-engine/tutorials/community/sentinel-2-s2cloudless

//   function maskS2sr(image) {
//     // Bits 10 and 11 are clouds and cirrus, respectively.
//     var cloudBitMask = ee.Number(2).pow(10).int();
//     var cirrusBitMask = ee.Number(2).pow(11).int();
//     // Get the pixel QA band.
//     var qa = image.select('QA60');
//     // All flags should be set to zero, indicating clear conditions.
//     var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
//         .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
//     // Return the masked image, scaled to TOA reflectance, without the QA bands.
//     return image.updateMask(mask)
//         .copyProperties(image, ["system:time_start"]);
//   }

//     var S2_cloudFilter_dry   = S2_SRdry.map(maskS2sr);
//     var S2_cloudFilter_flood = S2_SRflood.map(maskS2sr);
    
// // Get McFeeter's NDVI. Green-NIR/Green+NIR
// // < 0.3  = Non-Water
// // >= 0.3 = Water

//     var S2_clippeddry = S2_cloudFilter_dry.mosaic().clip(sindh);
//     var NDVIdry = S2_clippeddry.normalizedDifference(['B8','B4']).rename('NDVI');
//     // Rename band
//     var S2_cleandry = S2_clippeddry.addBands(NDVIdry);
    
//     var S2_clippedflood = S2_cloudFilter_flood.mosaic().clip(sindh);
//     var NDVIflood = S2_clippedflood.normalizedDifference(['B8','B4']).rename('NDVI');
//     // Rename band
//     var S2_cleanflood = S2_clippedflood.addBands(NDVIflood);

    
//     // S-2 Preview TODO: Fix bands, get NDVI and do reference samples
//     // False Colour Composite. NDVI - Red - Green
    
//   var vizParamsRGB = {
//   bands:['B4','B3','B2'],
//   max:3000
//   };
    
//     var vizParamsFCC = {
//   bands:['B8','B4','B3'],
//   max:3000
//   };
    
//   var vizParamsNDVI = {
//     bands:['NDVI'],
//     min: -1,
//     max: 1,
//     palette: ['Navy', 'Blue','BlueViolet', 'Magenta','Red','Orange','Yellow','Green','Lime']
//   };
//    // Map.addLayer(S2_cleandry, vizParamsRGB,'Sentinel-2 Before',false);
//     Map.addLayer(S2_cleandry, vizParamsFCC,'Sentinel-2 Before | FCC',false);
//    // Map.addLayer(S2_cleandry, vizParamsNDVI,'Sentinel-2 Before| NDVI',true);
    
//    // Map.addLayer(S2_cleanflood, vizParamsRGB,'Sentinel-2 Flood',false);
//     Map.addLayer(S2_cleanflood, vizParamsFCC,'Sentinel-2 Flood | FCC',false);
//     //Map.addLayer(S2_cleanflood, vizParamsNDVI,'Sentinel-2 Flood | NDVI',true);