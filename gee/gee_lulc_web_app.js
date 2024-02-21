/**************************************************************
Application of LU/LC changes in the Namibe Province, Angola. 
The app use as the main sources>
- Dynamic World (DW) dataset 
- Sentinel 2 Data

The application serve to the purpose of comparison land use and land cover changes from 2016 to 2023 based on Sentinel 2 data
Project Climate Resilience and Sustainable Water Management, PIN and FRESAN
Link to gee code: https://code.earthengine.google.com/4be09c87a3ab577a8e6f041282162add
******************************************************************/

//AOI
var namibeBoundary =  ee.Geometry.Polygon(
    [[[13.83562071063496,-17.4275303929221],
      [13.83562071063496,-13.39319000170853],
      [11.57243711688496,-13.39319000170853],
      [11.57243711688496,-17.4275303929221]]], null, false);


var geometry = namibeBoundary;

Map.centerObject(geometry, 10);

/*
* CONFIGURE THE IMAGERY
*/

var sentinelImages = {
'2016': createSentinelLayer('2023-01-01', '2023-10-30', geometry),
'2017': createSentinelLayer('2022-01-01', '2022-12-30', geometry),
'2018': createSentinelLayer('2021-01-01', '2021-12-30', geometry),
'2019': createSentinelLayer('2020-01-01', '2020-12-30', geometry),
'2020': createSentinelLayer('2019-01-01', '2019-12-30', geometry),
'2021': createSentinelLayer('2018-01-01', '2018-12-30', geometry),
'2022': createSentinelLayer('2017-01-01', '2017-12-30', geometry),
'2023': createSentinelLayer('2016-01-01', '2016-12-30', geometry),
};

var lulcImages = {
'2016': createLULCLayer('2023-01-01', '2023-10-30', geometry),
'2017': createLULCLayer('2022-01-01', '2022-12-30', geometry),
'2018': createLULCLayer('2021-01-01', '2021-12-30', geometry),
'2019': createLULCLayer('2020-01-01', '2020-12-30', geometry),
'2020': createLULCLayer('2019-01-01', '2019-12-30', geometry),
'2021': createLULCLayer('2018-01-01', '2018-12-30', geometry),
'2022': createLULCLayer('2017-01-01', '2017-12-30', geometry),
'2023': createLULCLayer('2016-01-01', '2016-12-30', geometry),
};

var probabilityImages = {
'2016': createProbabilityLayer('2023-01-01', '2023-10-30', geometry),
'2017': createProbabilityLayer('2022-01-01', '2022-12-30', geometry),
'2018': createProbabilityLayer('2021-01-01', '2021-12-30', geometry),
'2019': createProbabilityLayer('2020-01-01', '2020-12-30', geometry),
'2020': createProbabilityLayer('2019-01-01', '2019-12-30', geometry),
'2021': createProbabilityLayer('2018-01-01', '2018-12-30', geometry),
'2022': createProbabilityLayer('2017-01-01', '2017-12-30', geometry),
'2023': createProbabilityLayer('2016-01-01', '2016-12-30', geometry),
};

var ndviImages = {
'2016': calculateNDVI('2023-01-01', '2023-10-30', geometry),
'2017': calculateNDVI('2022-01-01', '2022-12-30', geometry),
'2018': calculateNDVI('2021-01-01', '2021-12-30', geometry),
'2019': calculateNDVI('2020-01-01', '2020-12-30', geometry),
'2020': calculateNDVI('2019-01-01', '2019-12-30', geometry),
'2021': calculateNDVI('2018-01-01', '2018-12-30', geometry),
'2022': calculateNDVI('2017-01-01', '2017-12-30', geometry),
'2023': calculateNDVI('2016-01-01', '2016-12-30', geometry),
};

var ndmiImages = {
'2016': calculateNDMI('2023-01-01', '2023-10-30', geometry),
'2017': calculateNDMI('2022-01-01', '2022-12-30', geometry),
'2018': calculateNDMI('2021-01-01', '2021-12-30', geometry),
'2019': calculateNDMI('2020-01-01', '2020-12-30', geometry),
'2020': calculateNDMI('2019-01-01', '2019-12-30', geometry),
'2021': calculateNDMI('2018-01-01', '2018-12-30', geometry),
'2022': calculateNDMI('2017-01-01', '2017-12-30', geometry),
'2023': calculateNDMI('2016-01-01', '2016-12-30', geometry),
};


/*
FUNCTIOS
- Function to mask clouds using the Sentinel-2 QA band.
- Construct a collection of corresponding Dynamic World and Sentinel-2 composites for inspection.
- Filter the DW and S2 collections by region and date.
*/

function maskS2clouds(image) {
var qa = image.select('QA60')

// Bits 10 and 11 are clouds and cirrus
var cloudBitMask = 1 << 10;
var cirrusBitMask = 1 << 11;

// Both flags should be set to zero, indicating clear conditions.
var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
           qa.bitwiseAnd(cirrusBitMask).eq(0))

// Return the masked and scaled data, without the QA bands.
return image.updateMask(mask).divide(10000)
    .select("B.*")
    .copyProperties(image, ["system:time_start"])
}

function createSentinelLayer(startDate, endDate, geometry) {
var s2 = ee.ImageCollection('COPERNICUS/S2_HARMONIZED') 
           .filterDate(startDate, endDate)
           .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 35))
           .filterBounds(geometry)
           .select(['QA60', 'B12', 'B8', 'B3','B4' ]) // Load only necessary bands
           .map(maskS2clouds);
           
return s2.median().clip(geometry).visualize({bands: ['B12', 'B8', 'B4'], min: 0, max: 0.35, gamma: 0.83});
}

function calculateNDVI(startDate, endDate, geometry) {
var sentinelData = ee.ImageCollection('COPERNICUS/S2_HARMONIZED') 
           .filterDate(startDate, endDate)
           .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 35))
           .filterBounds(geometry)
           .select(['QA60', 'B8', 'B4']) // Load only necessary bands
           .map(maskS2clouds);

// Calculate NDVI
var ndvi = sentinelData.map(function(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  return image.addBands(ndvi);
});

var ndviComposite = ndvi.median().clip(geometry).select('NDVI');

// Visualization parameters
var visParams = {
  min: -0.1,
  max: 0.7,
  palette: ['419bdf','ffffff', 'ce7e45',  'df923d', 'f1b555', 'fcd163', '99b718', '74a901',
  '66a000', '529400', '3e8601', '207401', '056201', '004c00', '023b01',
  '012e01', '011d01', '011301'] // pallete from modis ndvi, adjusted,'d8D8D8'
};

// Composite and clip
var composite = ndviComposite.clip(geometry);

// Visualize NDVI categories
return composite.visualize(visParams);
}

function calculateNDMI(startDate, endDate, geometry) {
var sentinelData = ee.ImageCollection('COPERNICUS/S2_HARMONIZED') 
           .filterDate(startDate, endDate)
           .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 35))
           .filterBounds(geometry)
           .select(['QA60', 'B8','B11']) // Load only necessary bands
           .map(maskS2clouds);

// Calculate NDMI
var ndmi = sentinelData.map(function(image) {
  var ndmi = image.normalizedDifference(['B8', 'B11']).rename('NDVI');
  return image.addBands(ndmi);
});

var ndmiComposite = ndmi.median().clip(geometry).select('NDVI');

// Visualization parameters
var visParams = {
  min: -0.24,
  max: 0.24,
  palette: ['red', 'yellow', 'green',  'blue'] 
};

// Composite and clip
var composite = ndmiComposite.clip(geometry);

// Visualize NDMI categories
return composite.visualize(visParams);
}

function createLULCLayer(startDate, endDate, geometry) {
var dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
           .filterDate(startDate, endDate)
           .filterBounds(geometry);
var classification = dw.select('label');
return classification.reduce(ee.Reducer.mode()).clip(geometry).visualize({min: 0, max: 8, palette: [
  '#419BDF', '#397D49', '#88B053', '#7A87C6', '#E49635', '#DFC35A',
  '#C4281B', '#A59B8F', '#B39FE1']});
}

function createProbabilityLayer(startDate, endDate, geometry) {
var dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
           .filterDate(startDate, endDate)
           .filterBounds(geometry);
var classification = dw.select('label');
var dwComposite = classification.reduce(ee.Reducer.mode()).clip(geometry);

var probabilityBands = [
  'water', 'trees', 'grass', 'flooded_vegetation', 'crops', 'shrub_and_scrub',
  'built', 'bare', 'snow_and_ice'
];

// Select probability bands
var probabilityCol = dw.select(probabilityBands); 

// Create an image with the average pixel-wise probability of each class across the time-period.
var meanProbability = probabilityCol.reduce(ee.Reducer.mean()); 
var projection = ee.Projection('EPSG:3857').atScale(10);
var meanProbability = meanProbability.setDefaultProjection(projection);

// Create the Top-1 Probability Hillshade.
var top1Probability = meanProbability.reduce(ee.Reducer.max());
var top1Confidence = top1Probability.multiply(100).int();
var hillshade = ee.Terrain.hillshade(top1Confidence).divide(255);

var dwVisParams = {
min: 0,
max: 8,
palette: [
  '#419BDF', '#397D49', '#88B053', '#7A87C6', '#E49635', '#DFC35A',
  '#C4281B', '#A59B8F', '#B39FE1'
  ]
};

// Colorize the classification image resulting hillshade
var rgbImage = dwComposite.visualize(dwVisParams).divide(255);
var probabilityHillshade = rgbImage.multiply(hillshade).clip(geometry);

return probabilityHillshade.visualize({min: 0, max: 0.8}); 
}

function updateMap(map, layerType, year) {
var layerImage;
switch (layerType) {
  case 'Sentinel 2':
    layerImage = sentinelImages[year];
    break;
  case 'NDVI':
    layerImage = ndviImages[year];
    break;
  case 'NDMI':
    layerImage = ndmiImages[year];
    break;
  case 'LU/LC Dynamic World':
    layerImage = lulcImages[year];
    break;
  case 'LU/LC Hilshade Probability':
    layerImage = probabilityImages[year];
    break;
  default:
    layerImage = null; //sentinelImages[year]; null;
}

if (layerImage) {
  var layer = ui.Map.Layer(layerImage);
  map.layers().set(0, layer);
}
}


/*
* SET UP THE MAPS AND SPLIT WIDGET
*/

// Create the left and right maps
var leftMap = ui.Map();
leftMap.setControlVisibility(false);
var rightMap = ui.Map();
rightMap.setControlVisibility(false);

leftMap.centerObject(geometry, 10);
rightMap.centerObject(geometry, 10);

// Define the addLayerSelector function
function addLayerSelector(mapToChange, defaultValue, position, layerSelectDropdown) {
var label = ui.Label('Choose a year');
var select = ui.Select({
  items: Object.keys(lulcImages), // Assuming this contains the correct years
  onChange: function(year) {
    var layerType = layerSelectDropdown.getValue();
    updateMap(mapToChange, layerType, year);
  }
});
select.setValue(Object.keys(lulcImages)[defaultValue], true);

var controlPanel = ui.Panel({widgets: [label, select], style: {position: position}});
mapToChange.add(controlPanel);

return select;
}


// Create a SplitPanel to hold the maps
var splitPanel = ui.SplitPanel({
firstPanel: leftMap,
secondPanel: rightMap,
orientation: 'horizontal',
wipe: true,
style: {stretch: 'both'}
});

// Add the SplitPanel to the UI root
ui.root.clear();
ui.root.add(splitPanel);


//DEFINE SIDE PANEL AND ITS CONTETS
// Create side panel and its elements
var sidePanel = ui.Panel({
style: { 
  width: '300px',
  position: 'top-left',
  padding: '8px'
}
});

// Add a header to the panel
var header = ui.Label({
value: 'LU/LC Analysis Tool',
style: { fontWeight: 'bold', fontSize: '26px', margin: '0 0 4px 0' }
});
sidePanel.add(header);

// Create a panel to hold the description and bullet points
var descriptionPanel = ui.Panel({
style: { 
  padding: '8px',
  margin: '4px 0'
},
layout: ui.Panel.Layout.Flow('vertical')
});
// Add the descriptive text
var descriptionText = 'Select layers and years to compare land use and land cover changes.'+ 
'Use the dropdown menus to configure the map display. The application offers:';
descriptionPanel.add(ui.Label({value: descriptionText, style: {fontSize: '13px', margin: '4px 0'}}));

// Add bullet points
var bulletPoints = [
'High-resolution satellite imagery analysis',
'Near-real-time LULC data visualization',
'Comparative analysis across different years'
];

bulletPoints.forEach(function(point) {
var bulletLabel = ui.Label({
  value: '• ' + point,
  style: {fontSize: '13px', margin: '2px 0', whiteSpace: 'pre'}
});
descriptionPanel.add(bulletLabel);
});

// Add the description panel to the side panel
sidePanel.add(descriptionPanel);


//ADD DROPDOWNS AND OTHER WIDGETS TO SIDEPANEL
// Dropdown for layer selection on left-side
var leftdescription = ui.Label({
value: 'Select layer for left-side',
style: { fontSize: '11px', margin: '10px 0 1px 0' }
});
sidePanel.add(leftdescription);
var leftLayerSelect = ui.Select({
items: ['Sentinel 2', 'NDVI', 'NDMI', 'LU/LC Dynamic World', 'LU/LC Hilshade Probability'],
value: 'Sentinel 2', // Set default value
placeholder: 'Select layer for left-side'
});
sidePanel.add(leftLayerSelect);

// Dropdown for layer selection on right-side
var rightdescription = ui.Label({
value: 'Select layer for right-side',
style: { fontSize: '11px', margin: '10px 0 1px 0' }
});
sidePanel.add(rightdescription);
var rightLayerSelect = ui.Select({
items: ['Sentinel 2', 'NDVI', 'NDMI', 'LU/LC Dynamic World', 'LU/LC Hilshade Probability'],
value: 'LU/LC Hilshade Probability', // Set default value
placeholder: 'Select layer for right-side'
});
sidePanel.add(rightLayerSelect);

// LEGEND
//----- start-----LEGEND FOR LU/LC
// Define the legend 
var legend = ui.Panel({
style: {
  padding: '8px 15px',
  position: 'bottom-center'
}
});

var legendTitle = ui.Label({
value: 'LU/LC Legend',
style: {fontWeight: 'bold', fontSize: '16px', margin: '0 0 4px 0'}
});
legend.add(legendTitle);

// Predefined list of labels and corresponding colors
var labels = ['water', 'trees', 'grass', 'crops',
  'shrub and scrub', 'built up and rocks', 'bare land']; // Add all labels
var colors = ['419bdf', '397d49', '88b053', 'e49635', 'dfc35a', 'c4281b',
  'a59b8f']; // Corresponding colors

labels.forEach(function(label, index) {
var colorBox = ui.Label({
  style: {
    backgroundColor: colors[index],
    // Ensure the color box is square
    padding: '8px',
    margin: '0 0 4px 0'
  }
});

var description = ui.Label({
  value: label,
  style: {margin: '0 0 4px 6px'}
});

legend.add(ui.Panel({
  widgets: [colorBox, description],
  layout: ui.Panel.Layout.Flow('horizontal')
}));
});

// Add the legend to the side panel
sidePanel.add(legend);
//----- end-----LEGEND FOR LU/LC

//----- start-----LEGEND FOR NDVI
// Create a panel for the NDVI legend
var ndviLegend = ui.Panel({
style: {
  padding: '8px 15px',
  position: 'bottom-center',
  stretch: 'horizontal'
},
layout: ui.Panel.Layout.Flow('horizontal')
});

// Add title to the NDVI legend
var ndviLegendTitle = ui.Label({
value: 'NDVI:',
style: {fontWeight: 'bold', fontSize: '16px', margin: '0 0 4px 0'}
});

ndviLegend.add(ndviLegendTitle);

// Define the colors used in the NDVI visualization
var ndviColors = ['419bdf', 'ce7e45', 'df923d', 'f1b555', 'fcd163', '99b718', 
                '74a901', '66a000', '529400', '3e8601', '207401', '056201', 
                '004c00', '023b01', '012e01', '011d01', '011301'];

// Function to create a color box
function createColorBox(color) {
return ui.Label({
  style: {
    backgroundColor: '#' + color,
    // Make the color box narrow and horizontal
    padding: '8px 2px',
    margin: '0 1px',
    fontSize: '0px' // Hide any text (e.g., label text)
  }
});
}

// Add color boxes to the legend
ndviColors.forEach(function(color) {
ndviLegend.add(createColorBox(color));
});

// Add NDVI legend to the side panel
sidePanel.add(ndviLegend);
//----- end-----LEGEND FOR NDVI

//----- start-----LEGEND FOR NDMI
// Create a panel for the NDMI legend
var ndmiLegend = ui.Panel({
style: {
  padding: '8px 15px',
  position: 'bottom-center',
  stretch: 'horizontal'
},
layout: ui.Panel.Layout.Flow('horizontal')
});

// Add title to the NDMI legend
var ndmiLegendTitle = ui.Label({  // Changed variable name to ndmiLegendTitle
value: 'NDMI:',
style: {fontWeight: 'bold', fontSize: '16px', margin: '0 0 4px 0'}
});

ndmiLegend.add(ndmiLegendTitle);  // Adding to ndmiLegend

// Define the colors used in the NDMI visualization
var ndmiColors = [
  'FF0000', // Red
  'FF3300',
  'FF6600',
  'FF9900', // Orange-Yellow
  'FFCC00',
  'FFFF00', // Yellow
  'CCFF00',
  '99FF00', // Yellow-Green
  '66FF00',
  '33FF00', // Green
  '00FF00',
  '00FF33',
  '00FF66', // Green-Cyan
  '00FF99',
  '00FFCC', // Cyan
  '00FFFF', // Cyan-Blue
  '0000FF'  // Blue
];

// Add color boxes to the NDMI legend
ndmiColors.forEach(function(color) {
ndmiLegend.add(createColorBox(color));
});

// Add NDMI legend to the side panel
sidePanel.add(ndmiLegend);
//----- end-----LEGEND FOR NDMI

//LAYER INFO
var layerDescriptions = {
'Sentinel 2': 'Sentinel-2 offers high-resolution, multi-spectral imaging for environmental monitoring, focusing on vegetation, soil, water, and coastal areas. This app utilizes its harmonized collection, ensuring data consistency for all scenes post January 25, 2022.',
'NDVI':'NDVI (Normalized Difference Vegetation Index) is a common calculation in remote sensing for assessing vegetation health and density.',
'NDMI':'The NDMI (Normalized Difference Moisture Index) uses NIR and SWIR bands to monitor changes '+
        'in leaf water content by analyzing moisture levels. The value range of the NDMI is -1 to 1. Negative values of NDMI (red) correspond to barren soil. '+
        'Values around zero (light green) generally correspond to water stress. High, positive values represent high canopy without water stress (blue).'+
        'caused by leaf internal structure and dry matter, offering a precise measurement of plant hydration.',
'LU/LC Dynamic World': 'Dynamic World provides a 10m real-time Land Use/Land Cover dataset with class probabilities for nine categories, using Sentinel-2 L1C data from June 27, 2015, to present. It emphasizes cloud-free observations with rigorous cloud and shadow masking.',
'LU/LC Hilshade Probability': 'The Hilshade Probability layer visualizes Dynamic World LU/LC probabilities using a hillshade effect, based on the highest per-pixel probability. It interprets class confidence as elevation, creating a unique landscape view of class distribution.'
};

var layerLinks = {
'Sentinel 2': 'https://sentinel.esa.int/web/sentinel/user-guides/sentinel-2-msi/resolutions/radiometric',
'NDVI':'https://www.sciencedirect.com/topics/earth-and-planetary-sciences/normalized-difference-vegetation-index',
'NDMI':'https://custom-scripts.sentinel-hub.com/sentinel-2/ndmi/',
'LU/LC Dynamic World': 'https://www.nature.com/articles/s41597-022-01307-4',
'LU/LC Hilshade Probability': 'https://www.nature.com/articles/s41597-022-01307-4'
};

// Function to update the layer information in the side panel
function updateLayerInfo(layerType) {
// Clear previous info
layerInfoPanel.clear();

// Add new description
var description = layerDescriptions[layerType];
layerInfoPanel.add(ui.Label({value: description, style: {fontSize: '13px', margin: '4px 0'}}));

// Add clickable link
var link = ui.Label('Read more', {color: 'blue', fontSize: '13px', margin: '4px 0'}, layerLinks[layerType]);
layerInfoPanel.add(link);
}

leftLayerSelect.onChange(function(layerType) {
var year = leftSelector.getValue();
updateMap(leftMap, layerType, year);
updateLayerInfo(layerType);
});

rightLayerSelect.onChange(function(layerType) {
var year = rightSelector.getValue();
updateMap(rightMap, layerType, year);
updateLayerInfo(layerType);
});
var layerInfoPanel = ui.Panel({
style: { 
  padding: '8px',
  margin: '10px 0'
}
});
sidePanel.add(layerInfoPanel);


// ADD THE SIDE PANEL TO THE UI
ui.root.insert(0, sidePanel);


/*
*LAYER SELECTOR
*/
// Create layer selectors for each map and link them to the dropdowns
var leftSelector = addLayerSelector(leftMap, 0, 'top-left', leftLayerSelect);
var rightSelector = addLayerSelector(rightMap, 0, 'top-right', rightLayerSelect);

// Set onChange for layer select dropdowns
leftLayerSelect.onChange(function(layerType) {
var year = leftSelector.getValue();
updateMap(leftMap, layerType, year);
});

rightLayerSelect.onChange(function(layerType) {
var year = rightSelector.getValue();
updateMap(rightMap, layerType, year);
});


// Define a default 
// Update both maps with the default layer 
updateMap(leftMap, 'Sentinel 2', '2016');
updateMap(rightMap, 'LU/LC Hilshade Probability', '2016');

// Update layer information for the default layer
updateLayerInfo('LU/LC Hilshade Probability');

// Link the maps together
var linker = ui.Map.Linker([leftMap, rightMap]);
leftMap.centerObject(geometry, 10);

