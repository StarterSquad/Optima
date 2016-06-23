define(['angular', 'jquery', './svg-to-png', 'underscore'], function (angular, $, svgToPng, _) {
  'use strict';

  return angular
      .module('app.mpld3-charts.export-helpers',[])
      .factory('exportHelpers', ['$http', function ($http) {

    var chartCssUrl = '/assets/css/chart.css';

    /**
     * Generate an image out of a SVG graph
     *
     * In order to achieve this a new SVG is created including styles.
     * This SVG element is used as data source inside an image which then
     * is used to draw the content on a canvas to save it as PNG.
     *
     * @param {element} el - DOM element containing the chart.
     * @param {function} callback - DOM element containing the chart.
     * @param {string} [type] - for 'blob' a PNG blob is returned, if undefined
     *                          a string of JPEG image data is returned
     */
    var generateGraphAsPngOrJpeg = function( el, callback, type ) {

      var elementId = $(el).find('[mpld3-chart]').attr('id');
      var $originalSvg = el.find('svg');
      var orginalWidth = $originalSvg.outerWidth();
      var orginalHeight = $originalSvg.outerHeight();
      var originalStyle = $originalSvg.attr('style') || '';
      var scalingFactor = 4.2;

      // in order to have styled graphs the css content used to render
      // graphs is retrieved & inject it into the svg as style tag
      var cssContentRequest = $http.get(chartCssUrl);
      cssContentRequest.success(function(chartStylesheetContent) {

        // It is needed to fetch all as mpld3 injects multiple style tags into the DOM
        var $styleTagContentList = $('style').map(function(index, style) {
          var styleContent = $(style).html();
          if (styleContent.indexOf('div#' + elementId) != -1) {
            return styleContent.replace(/div#/g, '#');
          } else {
            return styleContent;
          }
        });

        var styleContent = $styleTagContentList.get().join('\n');
        styleContent = styleContent + '\n' + chartStylesheetContent;

        // make sure we scale the padding and append it to the original styling
        // info: later declarations overwrite previous ones
        var svgInlineStyle = originalStyle + '; background:#fff; ' + svgToPng.scalePaddingStyle($originalSvg, scalingFactor);

        // create svg element
        var svg = svgToPng.createSvg(orginalWidth, orginalHeight, scalingFactor, svgInlineStyle, elementId);

        // add styles and content to the svg
        var styles = '<style>' + styleContent + '</style>';
        svg.innerHTML = styles + $originalSvg.html();

        // create img element with the svg as data source
        var svgXML = (new XMLSerializer()).serializeToString(svg);
        var tmpImage = document.createElement("img");
        tmpImage.width = orginalWidth * scalingFactor;
        tmpImage.height = orginalHeight * scalingFactor;
        tmpImage.src = "data:image/svg+xml;charset=utf-8,"+ svgXML;
        // console.log(svgXML);
        tmpImage.onload = function() {
          // draw image into canvas in order to convert it to a blob
          var canvas = document.createElement("canvas");
          canvas.width = orginalWidth * scalingFactor;
          canvas.height = orginalHeight * scalingFactor;
          var ctx = canvas.getContext("2d");
          ctx.drawImage(tmpImage, 0, 0);

          // Return the png either as a blob or data url
          if ( type == 'blob' ) {
            canvas.toBlob( callback );
          } else {
            var data = canvas.toDataURL('image/jpeg', 0.7);
            callback( data );
          }
        };
      }).error(function() {
        alert("Please releod and try again, something went wrong while generating the graph.");
      });
    };

    /**
     * Converts array of graph points to dictionary.
     */
    var graphToDictionary = function(graphPoints){
      var graphDict = {};
      _(graphPoints).each(function(point){
        graphDict[point[0]] = point[1];
      });
      return graphDict;
    };

    /**
     * Fills array according to reference dictionary keys,
     * Copies values from source dictionary if key found there, otherwise creates empty entries
     */
    var fillFromDictionary = function(reference, source){
      var result = [];
      _(reference).each(function(key){
        var nextEntry = (key in source)? source[key] : ''; // source[key] || '' would not catch '0' entries
        result.push(nextEntry);
      });
      return result;
    };

    /**
     * Returns the normalized data ready to export for a lineAreaScatter chart
     */
    var lineAreaScatterExport = function (graph){
      var exportable = {
        name: graph.options.title || 'Data',
        columns: []
      };

      // these fallback values are needed for cost-coverage since no legend should be
      // shown there, but they still need a title in the export
      // This certainly can be refactored to avoid the hardcoded column titles.
      var lineTitles = graph.options.legend? graph.options.legend : ["Line", "High", "Low"];

      // Collect and sort all the xPoints from the lines and the scatter data.
      // It's good enough to use the x data from the first line as all lines
      // have the same x data.
      var scatter = graphToDictionary(graph.data.scatter);
      var firstLine = graphToDictionary(_(graph.data.lines).first());
      var xOfPoints = {};
      xOfPoints.data = [];
      xOfPoints.title = graph.options.xAxis.axisLabel;
      xOfPoints.data.push.apply(xOfPoints.data, Object.keys(firstLine));
      xOfPoints.data.push.apply(xOfPoints.data, Object.keys(scatter));
      xOfPoints.data.sort();
      exportable.columns.push(xOfPoints);

      _(graph.data.lines).each(function(lineData, index) {

        var areaExists = graph.data.areas && graph.data.areas[index];
        // add low line
        if (areaExists) {
          var lowLine = graphToDictionary(graph.data.areas[index].lowLine);
          var yOfLowLinePoints = {};
          yOfLowLinePoints.title = 'Low ' + lineTitles[index];
          yOfLowLinePoints.data = fillFromDictionary(xOfPoints.data, lowLine);
          exportable.columns.push(yOfLowLinePoints);
        }

        // Collecting the Y of the points for the line
        var line = graphToDictionary(lineData);
        var yOfLinePoints = {};
        yOfLinePoints.title = lineTitles[index];
        yOfLinePoints.data = fillFromDictionary(xOfPoints.data, line);
        exportable.columns.push(yOfLinePoints);

        // add high line
        if (areaExists) {
          var highLine = graphToDictionary(graph.data.areas[index].highLine);
          var yOfHighLinePoints = {};
          yOfHighLinePoints.title = 'High ' + lineTitles[index];
          yOfHighLinePoints.data = fillFromDictionary(xOfPoints.data, highLine);
          exportable.columns.push(yOfHighLinePoints);
        }
      });

      var scatterPoints = {};
      scatterPoints.title = "scatter";
      scatterPoints.data = fillFromDictionary(xOfPoints.data, scatter);
      exportable.columns.push(scatterPoints);

      return exportable;
    };

    /**
     * Returns the normalized data ready to export for a stacked area chart
     */
    var areasExport = function (graph){
      var exportable = {
        name: graph.options.title || 'Data',
        columns: []
      };

      var areaTitles = graph.options.legend;

      // The X of the points are only sent in one column and we collect them from any of the lines
      var xOfPoints = {};
      xOfPoints.title = graph.options.xAxis.axisLabel;
      xOfPoints.data = _.map(graph.data.areas[0],function(point,j){ return point[0]; });
      exportable.columns.push(xOfPoints);

      _(graph.data.areas).each(function(areaData, index) {
        // Collecting the Y of the points for the area
        var yOfAreaPoints = {};
        yOfAreaPoints.title = areaTitles[index];
        yOfAreaPoints.data = _.map(areaData,function(point,j){ return point[1]; });
        exportable.columns.push(yOfAreaPoints);
      });

      return exportable;
    };

    /**
     * Returns the normalized data ready to export for a Radar Chart
     */
    var axesExport = function (graph){
      //x and y are not needed to be exported - they are just internal values to draw radar chart properly
      var exportable = {
        name: graph.options.title || 'Data',
        columns: []
      };

      var axisData = {};
      axisData.title = graph.radarAxesName;
      axisData.data = _.map(graph.data[0].axes, function(axis, j) { return axis.axis; });
      exportable.columns.push(axisData);

      _(graph.data).each(function(radarData, index) {
        var valueData = {};
        valueData.title = graph.options.legend[index];
        valueData.data = _.map(graph.data[index].axes, function(axis,j) { return axis.value; });
        exportable.columns.push(valueData);
      });

      return exportable;
    };

    /**
    * Returns the normalized data ready to export for a Radar Chart
    */
    var pieExport = function (graph){
      var exportable = {
        name: graph.options.title || 'Data',
        columns: []
      };

      _(graph.data.slices).each(function(slice, index) {
        var valueData = {};
        valueData.title = slice.label;
        valueData.data = [slice.value];
        exportable.columns.push(valueData);
      });

      return exportable;
    };

    /**
    * Returns the normalized data ready to export for a stacked bar charts
    */
    var stackedBarsExport = function (chart){
      var exportable = {
        name: chart.options.title,
        columns: []
      };

      var xData = {};
      xData.title = chart.options.xAxis.axisLabel;
      xData.data = _(chart.data.bars).map(function(bar) { return bar[0]; });
      exportable.columns.push(xData);

      _(chart.options.legend).each(function(title, index) {
        var yData = {};
        yData.title = title;

        yData.data = _(chart.data.bars).map(function(bar) { return bar[1][index]; });
        exportable.columns.push(yData);
      });

      return exportable;
    };

    /**
     * Returns the normalized data ready to export from any type of chart.
     *
     * @param {object} chart - an object which must contain the chart data.
     */
    var getExportableFrom = function (chart){
      if(!chart.data) { return null; }

      if(_.isEqual(Object.keys(chart.data),["lines", "scatter", "areas"]) ||
         _.isEqual(Object.keys(chart.data),["lines", "scatter", "limits"]) ||
         _.isEqual(Object.keys(chart.data),["lines", "scatter"])) {
        return lineAreaScatterExport(chart);
      }
      if(_.isEqual(Object.keys(chart.data),["areas"])) { return areasExport(chart); }
      if(_.isEqual(Object.keys(chart.data),["slices"])) { return pieExport(chart); }
      if(_.isEqual(Object.keys(chart.data),["bars"])) { return stackedBarsExport(chart); }
      if(_.isEqual(chart.data[0] && Object.keys(chart.data[0]),["axes"])) { return axesExport(chart); }
      return null;
    };

    /**
     * Returns the normalized data of an Mpld3 chart ready to export.
     *
     * @param {object} chart - an object which must contain the Mpld3 chart data.
     */
    var getMpld3ExportableFrom = function (chart) {
      if(!chart.data) { return null; }

      var exportable = {
        name: 'Data',
        columns: []
      };

      var xColumnData = [];
      _(chart.data).each(function(entry) {
        _(entry).each(function(point) {
          xColumnData.push(point[0]);
        });
      });
      xColumnData.sort();
      exportable.columns.push({
        title: chart.axes[0].texts[0].text,
        data: xColumnData
      });

      _(chart.data).each(function(entry) {
        var column = {
          title: chart.axes[0].texts[1].text
        };
        column.data = fillFromDictionary(xColumnData, _.object(entry));
        exportable.columns.push(column);
      });

      return exportable;
    };

    return {
      generateGraphAsPngOrJpeg: generateGraphAsPngOrJpeg,
      getExportableFrom: getExportableFrom,
      getMpld3ExportableFrom: getMpld3ExportableFrom
    };
  }]);
});
