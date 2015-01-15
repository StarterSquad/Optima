define([
    'angular',
    'ui.router',
    '../../config',
    '../resources/model',
    '../ui/type-selector/index',
    '../common/export-all-charts',
    '../common/export-all-data',
    '../common/graph-type-factory',
    '../common/update-checkbox-on-click-directive'
], function (angular) {
    'use strict';

    return angular.module('app.analysis', [
        'app.constants',
        'app.export-all-charts',
        'app.export-all-data',
        'app.resources.model',
        'app.ui.type-selector',
        'ui.router',
        'app.common.graph-type',
        'app.common.update-checkbox-on-click'
    ]).config(function ($stateProvider) {
        $stateProvider
            .state('analysis', {
                url: '/analysis',
                abstract: true,
                template: '<div ui-view></div>'
            })
            .state('analysis.scenarios', {
                url: '/scenarios',
                templateUrl: 'js/modules/analysis/scenarios.html' ,
                controller: 'AnalysisScenariosController',
                resolve: {
                  info: function($http) {
                    return $http.get('/api/project/info');
                  },
                  meta: function (Model) {
                    return Model.getParametersDataMeta().$promise;
                  },
                  scenarioParamsResponse: function($http, info) {
                    return $http.get('/api/analysis/scenarios/params');
                  },
                  scenariosResponse: function($http, info) {
                    return $http.get('/api/analysis/scenarios/list');
                  }
                }
            })
            .state('analysis.optimization', {
                url: '/optimization',
                templateUrl: 'js/modules/analysis/optimization.html' ,
                controller: 'AnalysisOptimizationController',
                resolve: {
                  meta: function (Model) {
                    return Model.getParametersDataMeta().$promise;
                  }
                }
            });
    });

});
