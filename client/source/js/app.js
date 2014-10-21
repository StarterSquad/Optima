define([
  'angular',
  'angular-nvd3',
  'ui.bootstrap',
  'ui.router',
  './config',
  './modules/about/index',
  './modules/analyses/index',
  './modules/common/save-graph-as-directive',
  './modules/d3-charts/index',
  './modules/graphs/index',
  './modules/help/index',
  './modules/home/index',
  './modules/model/index',
  './modules/optimization/index',
  './modules/playground/index',
  './modules/project/index',
  './modules/import-export/index',
  './modules/results/index',
  './modules/ui/index'
], function (angular) {
  'use strict';

  return angular.module('app', [
    'app.about',
    'app.analyses',
    'app.constants',
    'app.d3-charts',
    'app.graphs',
    'app.help',
    'app.home',
    'app.import-export',
    'app.model',
    'app.optimization',
    'app.playground',
    'app.project',
    'app.results',
    'app.save-graph-as',
    'app.ui',
    'nvd3',
    'ui.bootstrap',
    'ui.router'
  ]).config(function ($urlRouterProvider) {
    $urlRouterProvider.otherwise('/');
  });

});
