/**
 * loads sub modules and wraps them up into the main module.
 * This should be used for top-level module definitions only.
 */
define([
  'angular',
  'ui.router',
  './config',
  './modules/graphs/index',
  './modules/home/index',
  './modules/ui/index'
], function (angular) {
  'use strict';

  return angular.module('app', [
    'app.constants',
    'app.graphs',
    'app.home',
    'app.ui',
    'ui.router'
  ]).config(['$urlRouterProvider', function ($urlRouterProvider) {
    $urlRouterProvider.otherwise('/');
  }]);

});
