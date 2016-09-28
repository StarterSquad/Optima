/**
 * app.home defines the module and controller for the home URI of the app
 * Makes the app to work either on a project that might be set or a defalut welcome view.
 */
define([
  'angular',
  'ui.router',
  '../common/active-project-service',
  '../../config'
], function (angular) {
  'use strict';

  return angular.module('app.home', [
    'app.constants',
    'app.active-project',
    'ui.router'
  ]).config(function ($stateProvider) {
    $stateProvider
      .state('home', {
        url: '/',
        templateUrl: 'js/modules/home/home.html',
        controller: 'ProjectOpenController',
        resolve: {
          projects: function (projectApiService) {
            return projectApiService.getProjectList();
          },
        }
      });
      // .state('home', {
      //   url: '/',
      //   templateUrl: './js/modules/home/home.html',
      //   controller: 'HomeController',
      //   resolve: {
      //     project: function (projectApiService) {
      //       return projectApiService.getActiveProject();
      //     }
      //   }
      // });
  });

});
