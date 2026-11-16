/*******************************************************************
 * ASPRS GeoChallenge 2026 - Pine Gulch Fire, Colorado
 * Download pre- and post-fire Sentinel-2 L2A images (12 bands)
 * Extent: MTBS burn-area shapefile, event CO3933610852620200731
 *******************************************************************/

// ================= 1. SETTINGS =================
var FIRE_ASSET = 'projects/saz001/assets/pine_gulch_mtbs';
var SITE       = 'pinegulch';
var BUFFER_M   = 2000;                 // margin around the perimeter, metres
var PRE_START  = '2020-06-25';
var PRE_END    = '2020-07-31';         // ignition day; end date is exclusive
var POST_START = '2020-09-25';
var POST_END   = '2020-10-25';
var MAX_CLOUD  = 40;                   // scene-level cloud filter, percent
var CRS        = 'EPSG:32612';         // WGS 84 / UTM zone 12N
var SCALE_M    = 20;                   // output pixel size, metres
var DRIVE_DIR  = 'GeoChallenge';

// ================= 2. AREA OF INTEREST =================
var fire = ee.FeatureCollection(FIRE_ASSET);
var aoi  = fire.geometry().bounds().buffer(BUFFER_M).bounds();
print('Features in shapefile (expect 1):', fire.size());
print('Bounding box:', aoi);

// ================= 3. SENTINEL-2 COMPOSITES =================
var BANDS = ['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11','B12'];

// Keep only clear pixels using the scene classification layer (SCL)
function maskS2(img) {
  var scl = img.select('SCL');
  var clear = scl.neq(0)     // no data
    .and(scl.neq(1))         // saturated or defective
    .and(scl.neq(3))         // cloud shadow
    .and(scl.neq(8))         // cloud, medium probability
    .and(scl.neq(9))         // cloud, high probability
    .and(scl.neq(10))        // thin cirrus
    .and(scl.neq(11));       // snow or ice
  return img.updateMask(clear).select(BANDS);
}

// Median of all clear observations in a date window
function composite(start, end, label) {
  var col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', MAX_CLOUD));
  print(label + ' scenes (' + start + ' to ' + end + '):', col.size());
  return col.map(maskS2).median().clip(aoi);
}

var pre  = composite(PRE_START,  PRE_END,  'Pre-fire');
var post = composite(POST_START, POST_END, 'Post-fire');

// ================= 4. VISUAL CHECK =================
var swirVis = {bands: ['B12','B8A','B4'], min: 0, max: 4000};
Map.centerObject(aoi, 10);
Map.addLayer(pre,  swirVis, 'Pre-fire (SWIR-NIR-Red)');
Map.addLayer(post, swirVis, 'Post-fire (SWIR-NIR-Red)');
Map.addLayer(fire.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'MTBS perimeter');
Map.addLayer(ee.Image().paint(aoi, 1, 2), {palette: ['yellow']}, 'Bounding box');

// ================= 5. EXPORT TO GOOGLE DRIVE =================
function exportImage(img, label) {
  Export.image.toDrive({
    image: img.toUint16(),
    description: SITE + '_S2_' + label,
    folder: DRIVE_DIR,
    region: aoi,
    scale: SCALE_M,
    crs: CRS,
    maxPixels: 1e10
  });
}
exportImage(pre,  'pre');
exportImage(post, 'post');