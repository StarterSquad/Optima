define(['./module', './scale-helpers', 'angular', './stacked-bar-chart'], function (module, scaleHelpers, angular, StackedBarChart) {
  'use strict';

  module.directive('stackedBarChart', function (d3Charts) {
    var svg;

    var defaultColors = [ '__light-blue', '__blue', '__violet', '__green',
    '__light-green', '__gray', '__red' ];

    /**
     * Draw the stacked bar chart
     */
    var drawGraph = function (data, options, rootElement) {
      options = d3Charts.adaptOptions(options);

      // to prevent creating multiple graphs we want to remove the existing svg
      // element before drawing a new one.
      if (svg) {
        rootElement.find("svg").remove();
      }

      var dimensions = {
        height: options.height,
        width: options.width
      };

      svg = d3Charts.createSvg(rootElement[0], dimensions, options.margin);

      var chartSize = {
        width: options.width - options.margin.left - options.margin.right,
        height: options.height - options.margin.top - options.margin.bottom
      };

      // Define svg groups
      var chartGroup = svg.append('g').attr('class', 'chart_group');
      var axesGroup = svg.append('g').attr('class', 'axes_group');
      var headerGroup = svg.append('g').attr('class', 'header_group');

      var chart = new StackedBarChart(chartGroup, chartSize, data.bars, options.linesStyle);
      chart.draw();
      d3Charts.drawTitleAndLegend(svg, options, headerGroup);

      options.yAxis.tickFormat = function (value) {
        var format = scaleHelpers.evaluateTickFormat(0, chart.yMax());
        return scaleHelpers.customTickFormat(value, format);
      };

      d3Charts.drawAxes(
        chart.scales(),
        options,
        axesGroup,
        chartSize
      );

    };

    return {
      scope: {
        data: '=',
        options: '='
      },
      link: function (scope, element) {
        // before this change, all the graphs were redrawn three times
        scope.$watchCollection('[data,options]', function() {
          drawGraph(scope.data, angular.copy(scope.options), element);
        });
      }
    };
  });
});
