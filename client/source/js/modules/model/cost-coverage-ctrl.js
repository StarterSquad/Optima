define(['./module', 'underscore'], function (module, _) {
  'use strict';

  module.controller('ModelCostCoverageController', function ($scope, $http,
    $state, meta, info, modalService, programs) {

    var plotTypes, effectNames;

    var initialize =function () {
      $scope.meta = meta;
      $scope.chartsForDataExport = [];
      $scope.titlesForChartsExport = [];

      // show message "calibrate the model" and disable the form elements
      $scope.projectInfo = info;
      $scope.needData = !$scope.projectInfo.has_data;
      $scope.cannotCalibrate = !$scope.projectInfo.can_calibrate;
      $scope.notReady = $scope.needData || $scope.cannotCalibrate;

      $scope.optionsErrorMessage = 'To define a cost-coverage curve, values must be provided in the first three text boxes.';
      $scope.all_programs = programs;

      if ( !$scope.needData ) {
        $scope.initializePrograms();
        $scope.selectedProgram = $scope.programs[0];
        $scope.displayedProgram = null;

        $scope.coParams = [];

        $scope.hasCostCoverResponse = false;
      }

      // model parameters
      $scope.defaultSaturationCoverageLevel = 90;
      $scope.defaultKnownCoverageLevel = 60;
      $scope.defaultKnownFundingValue = 400000;
      $scope.defaultScaleUpParameter = 1;
      $scope.defaultNonHivDalys = 0;
      $scope.defaultXAxisMaximum = 1000000;
      $scope.behaviorWithoutMin = 0.3;
      $scope.behaviorWithoutMax = 0.5;
      $scope.behaviorWithMin = 0.7;
      $scope.behaviorWithMax = 0.9;
      $scope.xAxisMaximum = undefined;
      $scope.saturationCoverageLevel = undefined;
      $scope.knownCoverageLevel = undefined;
      $scope.knownFundingValue = undefined;
      $scope.scaleUpParameter = undefined;
      $scope.nonHivDalys = undefined;
      $scope.displayCost = 1;

      plotTypes = ['plotdata', 'plotdata_cc', 'plotdata_co'];

      resetGraphs();
    };

    /**
     * Redirects the user to View & Calibrate screen.
     */
    $scope.gotoViewCalibrate = function() {
      $state.go('model.view');
    };

    /**
    * Creates the models of the programs for this controller.
    * If the backend do not present values for the categories, we'll use 'Others' as default.
    */
    $scope.initializePrograms = function () {
      $scope.programs =  _(meta.progs.long).map(function (name, index) {
        var acronym = meta.progs.short[index];
        return {
          name: name,
          acronym: acronym,
          category: 'Other', // it will be read from project_info, once it is synced with meta.programs
          ccparams: $scope.all_programs[acronym].ccparams,
          ccplot: $scope.all_programs[acronym].ccplot
        };
      });
      /** Dec 26 2014
       * fix/306-2-fix-plotting-of-default-ccocs
       * Default null value for selectedProgram
       */
      $scope.programs.unshift({name:'-- No program selected --',category:null, acronym:null});
    };


    var resetGraphs= function () {
      $scope.graphs = {
        plotdata: [],
        plotdata_cc: {},
        plotdata_co: []
      };
    };

    var getLineScatterOptions = function (options, xLabel, yLabel) {
      var defaults = {
        width: 300,
        height: 200,
        margin: {
          top: 20,
          right: 15,
          bottom: 40,
          left: 60
        },
        xAxis: {
          axisLabel: xLabel || 'X'
        },
        yAxis: {
          axisLabel: yLabel || 'Y'
        }
      };

      return _(angular.copy(defaults)).extend(options);
    };

    /* Methods
     ========= */

    /**
     * Calculates graphs objects of types plotdata and plotdata_co
     * returns ready to draw Graph object
     * @param graphData - api reply
     * @returns {{options, data: {lines: Array, scatter: Array}}}
     */
    var setUpPlotdataGraph = function (graphData) {

      var graph = {
        options: getLineScatterOptions({
          width: 300,
          height: 200,
          margin: {
            top: 20,
            right: 15,
            bottom: 40,
            left: 60
          },
          linesStyle: ['__blue', '__black __dashed', '__black __dashed']
        }, graphData.xlabel, graphData.ylabel),
        data: {
          lines: [],
          scatter: []
        },
        title: graphData.title
      };

      // quit if data is empty - empty graph placeholder will be displayed
      if (graphData.ylinedata) {

        var numOfLines = graphData.ylinedata.length;

        _(graphData.xlinedata).each(function (x, index) {
          var y = graphData.ylinedata;
          for (var i = 0; i < numOfLines; i++) {
            if (!graph.data.lines[i]) {
              graph.data.lines[i] = [];
            }

            graph.data.lines[i].push([x, y[i][index]]);
          }
        });
      }

      _(graphData.xscatterdata).each(function (x, index) {
        var y = graphData.yscatterdata;

        if (y[index]) {
          graph.data.scatter.push([x, y[index]]);
        }
      });

      // set up the data limits
      graph.data.limits = [
        [graphData.xlowerlim, graphData.ylowerlim],
        [graphData.xupperlim, graphData.yupperlim]
      ];


      return graph;
    };

    /**
     * Generates ready to plot graph for a cost coverage.
     */
    var prepareCostCoverageGraph = function (data) {
      var graph = {
        options: getLineScatterOptions({}, data.xlabel, data.ylabel),
        data: {
          // there is a single line for that type
          lines: [[]],
          scatter: []
        }
      };

      _(data.xlinedata).each(function (x, index) {
        var y = data.ylinedata;
        graph.data.lines[0].push([x, y[index]]);
      });

      _(data.xscatterdata).each(function (x, index) {
        var y = data.yscatterdata;

        if (y[index]) {
          graph.data.scatter.push([x, y[index]]);
        }
      });

      // set up the data limits
      graph.data.limits = [
        [data.xlowerlim, data.ylowerlim],
        [data.xupperlim, data.yupperlim]
      ];
      return graph;
    };

    /**
     * Receives graphs data with plot type to calculate, calculates all graphs
     * of given type and writes them to $scope.graphs[type] except for the
     * cost coverage graph which will be written to $scope.ccGraph
     *
     * @param data - usually api request with graphs data
     * @param type - string
     */
    var prepareGraphsOfType = function (data, type) {
      if (type === 'plotdata_cc') {
        $scope.ccGraph = prepareCostCoverageGraph(data);
        $scope.ccGraph.title = $scope.displayedProgram.name;
      } else if (type === 'plotdata' || type === 'plotdata_co') {
        _(data).each(function (graphData) {
          $scope.graphs[type].push(setUpPlotdataGraph(graphData));
        });
      }
    };

    var setUpCOParamsFromEffects = function (effectNames) {
      $scope.coParams = _(effectNames).map(function (effect) {
        return [
          (effect[2] && effect[2][0])? effect[2][0] : null,
          (effect[2] && effect[2][1])? effect[2][1] : null,
          (effect[2] && effect[2][2])? effect[2][2] : null,
          (effect[2] && effect[2][3])? effect[2][3] : null
        ];
      });
    };

    $scope.convertFromPercent = function (value) {
      if (typeof value !== "number" || isNaN(value)) {
        return NaN;
      }
      return value / 100;
    };

    $scope.costCoverageParams = function () {
      return [
        $scope.convertFromPercent($scope.saturationCoverageLevel),
        $scope.convertFromPercent($scope.knownCoverageLevel),
        $scope.knownFundingValue,
        $scope.scaleUpParameter,
        $scope.nonHivDalys
      ];
    };

    var ccPlotParams = function() {
      if ($scope.xAxisMaximum) {
        var years = [];
        if ($scope.displayCost == 2 && $scope.displayYear) {
          years = [1, [parseInt($scope.displayYear, 10)]];
        } else {
          years = [0, []];
        }
        return [$scope.xAxisMaximum, years];
      } else {
        return [];
      }
    };

    /**
     * Returns the current parameterised plot model.
     */
    var getPlotModel = function() {
      return {
        progname: $scope.selectedProgram.acronym,
        ccparams: $scope.costCoverageParams(),
        coparams: [],
        ccplot: ccPlotParams()
      };
    };

    /**
     * Returns true if all of the elements in an array are defined or not null
     */
    var hasAllElements = function(params) {
      return params && params.length && _(params).every(function(item) { return item; });
    };

    /**
     * Returns true if all of the elements in an array are undefined, null or NaN
     */
    var hasOnlyInvalidEntries = function(params) {
      return params.every(function(item) {
        return item === undefined || item === null || typeof item === "number" && isNaN(item);
      });
    };

    $scope.areValidParams = function (params) {
      return hasAllElements(params) || hasOnlyInvalidEntries(params);
    };

    var areCCParamsValid = function (params) {
      return $scope.areValidParams(params.slice(0, 3));
    };

    $scope.hasValidCCParams = function() {
      return !$scope.hasCostCoverResponse || areCCParamsValid($scope.costCoverageParams());
    };

    /**
     * Update current program ccparams based on the selected program.
     *
     * This function is supposed to be called before Draw / Redraw / Save.
     */
    var updateCCParams = function(model) {
      if (model.ccparams) {
        $scope.selectedProgram.ccparams = model.ccparams;
        $scope.all_programs[$scope.selectedProgram.acronym].ccparams = model.ccparams;
      }
      if (model.ccplot) {
        $scope.selectedProgram.ccplot = model.ccplot;
      }
    };

    /**
     * Retrieve and update graphs based on the provided plot models.
     */
    var retrieveAndUpdateGraphs = function (model) {
      // validation on Cost-coverage curve plotting options
      if ( !areCCParamsValid(model.ccparams) || $scope.ccForm.valid){
        modalService.inform(
          function () {},
          'Okay',
          $scope.optionsErrorMessage,
          'Error!'
        );
        return;
      }

      // stop further execution and return in case of null selectedProgram
      if ( $scope.selectedProgram.acronym === null ) {
        return;
      }

      // clean up model by removing unnecessary parameters
      if (_.isEmpty(model.ccparams) || hasOnlyInvalidEntries(model.ccparams.slice(0,3))) {
        delete model.ccparams;
      }

      if (_.isEmpty(model.coparams) || hasOnlyInvalidEntries(model.coparams)) {
        delete model.coparams;
      }

      // update current program ccparams,if applicable
      updateCCParams(model);

      $http.post('/api/model/costcoverage', model).success(function (response) {
        if (response.status === 'OK') {

          $scope.displayedProgram = angular.copy($scope.selectedProgram);
          effectNames = response.effectnames;
          setUpCOParamsFromEffects(response.effectnames);
          $scope.hasCostCoverResponse = true;

          resetGraphs();
          _(plotTypes).each(function (plotType) {
            prepareGraphsOfType(response[plotType], plotType);
          });
        }
      });
    };

    $scope.changeProgram = function() {
      if($scope.hasCostCoverResponse === true) {
        $scope.hasCostCoverResponse = false;
      }
      if (hasAllElements($scope.selectedProgram.ccparams.slice(0,3))) {
        $scope.saturationCoverageLevel = $scope.selectedProgram.ccparams[0]*100;
        $scope.knownCoverageLevel = $scope.selectedProgram.ccparams[1]*100;
        $scope.knownFundingValue = $scope.selectedProgram.ccparams[2];
        $scope.scaleUpParameter = $scope.selectedProgram.ccparams[3];
        $scope.nonHivDalys = $scope.selectedProgram.ccparams[4];
      } else {
        $scope.saturationCoverageLevel = undefined;
        $scope.knownCoverageLevel = undefined;
        $scope.knownFundingValue = undefined;
        $scope.scaleUpParameter = undefined;
        $scope.nonHivDalys = undefined;
      }
      if ($scope.selectedProgram.ccplot && $scope.selectedProgram.ccplot.length==2) {
        $scope.xAxisMaximum = $scope.selectedProgram.ccplot[0];
        var years = $scope.selectedProgram.ccplot[1][1];
        if (years.length > 0) {
          $scope.displayYear = years[0];
          $scope.displayCost = 2;
        } else {
          $scope.displayCost = 1;
          $scope.displayYear = undefined;
        }
      } else {
        $scope.displayCost = 1;
        $scope.displayYear = undefined;
        $scope.xAxisMaximum = undefined;
      }

      $scope.generateCurves();
    };

    /**
     * Retrieve and update graphs based on the current plot models.
     */
    $scope.generateCurves = function () {
      var model = getPlotModel();
      if ($scope.hasCostCoverResponse) {
        model.all_coparams = $scope.coParams;
        model.all_effects = effectNames;
      }
      retrieveAndUpdateGraphs(model);
    };

    $scope.uploadDefault = function () {
      var message = 'Upload default cost-coverage-outcome curves will be available in a future version of Optima. We are working hard in make it happen for you!';
      modalService.inform(
        function () {},
        'Okay',
        message,
        'Thanks for your interest!'
      );
    };

    /**
     * Retrieve and update graphs based on the current plot models.
     *
     * The plot model gets saved in the backend.
     */
    $scope.saveModel = function () {
      var model = getPlotModel(model);
      model.doSave = true;
      model.all_coparams = $scope.coParams;
      model.all_effects = effectNames;
      retrieveAndUpdateGraphs(model);
    };

    /**
     * Retrieve and update graphs based on the current plot models.
     *
     * The plot model gets reverted in the backend.
     */
    $scope.revertModel = function () {
      var model = getPlotModel(model);
      model.doRevert = true;
      retrieveAndUpdateGraphs(model);
    };

    /**
     * POST /api/model/costcoverage/effect
     *   {
     *     "progname":<chosen progname>
     *     "effect":<effectname for the given row>,
     *     "ccparams":<ccparams>,
     *     "coparams":<coprams from the corresponding coparams block>
     *   }
     */
    $scope.updateCurve = function (graphIndex) {
      var model = getPlotModel();
      model.coparams = $scope.coParams[graphIndex];
      model.effect =  effectNames[graphIndex];
      if ( !$scope.areValidParams(model.coparams) ){
        modalService.inform(
          function () {},
          'Okay',
          $scope.optionsErrorMessage,
          'Error!'
        );
        return;
      }

      // clean up model by removing unnecessary parameters
      if (_.isEmpty(model.ccparams) || hasOnlyInvalidEntries(model.ccparams)) {
        delete model.ccparams;
      }

      if (_.isEmpty(model.coparams) || hasOnlyInvalidEntries(model.coparams)) {
        delete model.coparams;
      }

      // update current program ccparams, if applicable
      updateCCParams(model);

      $http.post('/api/model/costcoverage/effect', model).success(function (response) {
        $scope.graphs.plotdata[graphIndex] = setUpPlotdataGraph(response.plotdata);
        $scope.graphs.plotdata_co[graphIndex] = setUpPlotdataGraph(response.plotdata_co);
        effectNames[graphIndex]=response.effect;
      });
    };

    /**
     * Collects all existing charts in the $scope.chartsForDataExport variable.
     * In addition all titles are gatherd into titlesForChartsExport. This is
     * needed since the cost coverage graphs have no title on the graphs.
     */
    var updateDataForExport = function() {
      $scope.chartsForDataExport = [];
      $scope.titlesForChartsExport = [];

      if ( $scope.ccGraph) {
        $scope.chartsForDataExport.push($scope.ccGraph);
        $scope.titlesForChartsExport.push($scope.ccGraph.title);
      }

      var charts = _(_.zip($scope.graphs.plotdata, $scope.graphs.plotdata_co)).flatten();
      _( charts ).each(function (chart,index) {
        $scope.chartsForDataExport.push(chart);
        $scope.titlesForChartsExport.push(chart.title);
      });
    };

    $scope.$watch('graphs', updateDataForExport, true);
    $scope.$watch('ccGraph', updateDataForExport, true);

    /**
     * Retrieve and update graphs based on the current plot models only if the graphs are already rendered
     * by pressing the draw button.
     */
    $scope.updateCurves =  _.debounce(function() { // debounce a bit so we don't update immediately
      if($scope.CostCoverageForm.$valid && $scope.hasCostCoverResponse === true) {
       $scope.generateCurves();
      }
    }, 500);

    initialize();

  });
});
