define([
  './module'
], function (module) {
  'use strict';
  module.controller('GraphsController', function ($scope, dataMocks, CONFIG) {

    // GRAPH 5
    // =======
    $scope.options5 = {
      height: 300,
      width: 900,
      margin: {
        top: 20,
        right: 20,
        bottom: 60,
        left: 60
      },
      xAxis: {
        axisLabel: 'Axis X',
        tickFormat: function (d) {
          return d3.format('d')(d);
        }
      },
      yAxis: {
        axisLabel: 'Axis Y',
        tickFormat: function (d) {
          return d3.format(',.2f')(d);
        }
      }
    };

    dataMocks.lineScatterError().$promise.then(function (data) {
      $scope.data5 = data;
    });



    // GRAPH 6
    // =======

    $scope.options6 = {
      height: 300,
      width: 900,
      margin: {
        top: 20,
        right: 20,
        bottom: 60,
        left: 60
      },
      xAxis: {
        axisLabel: 'Axis X',
        tickFormat: function (d) {
          return d3.format('d')(d);
        }
      },
      yAxis: {
        axisLabel: 'Axis Y',
        tickFormat: function (d) {
          return d3.format(',.2f')(d);
        }
      }
    };

    dataMocks.lineScatterArea().$promise.then(function (data) {
      $scope.data6 = data;
    });

    // GRAPH 7
    // =======

    $scope.options7 = {
      height: 300,
      width: 900,
      margin: {
        top: 20,
        right: 20,
        bottom: 60,
        left: 60
      },
      xAxis: {
        axisLabel: 'Axis X',
        tickFormat: function (d) {
          return d3.format('d')(d);
        }
      },
      yAxis: {
        axisLabel: 'Axis Y',
        tickFormat: function (d) {
          return d3.format(',.2f')(d);
        }
      }
    };

    $scope.data7 = {
      'areas': [
        [[2001, 1], [2002, 3], [2004, 5], [2005, 6], [2006, 7]],
        [[2001, 2], [2002, 1], [2004, 2], [2005, 2], [2006, 4]],
        [[2001, 1], [2002, 2], [2004, 20], [2005, 2], [2006, 1]]
      ]
    };

    // GRAPH 8
    // =======

    $scope.options8 = {
      height: 300,
      width: 300,
      margin: {
        top: 20,
        right: 20,
        bottom: 60,
        left: 60
      }
    };

    $scope.data8 = {
      slices: [
        {label: '<5', value: 2704659},
        {label: '5-13', value: 4499890},
        {label: '14-17', value: 4499890},
        {label: 'Hello my long text A', value: 5},
        {label: 'Ohhh, another long label', value: 10}
      ]
    };

    // Stacked bar chart
    // =======

    var linesStyle = [ '__light-blue', '__blue', '__violet', '__green',
    '__light-green', '__gray', '__red' ];

    $scope.options9 = {
      linesStyle: linesStyle,
      legend: ['AA', 'BB', 'CC'],
      height: 300,
      width: 300,
      margin: {
        top: 20,
        right: 20,
        bottom: 60,
        left: 60
      }
    };

    $scope.data9 = {
      bars: [
       [2001, [552339, 259034, 450818]],
       [2002, [85640, 42153, 74257]],
       [2003, [592339, 559034, 250818]]
      ]
    };

  });
});
