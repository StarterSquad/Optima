define([
  'angular',
  'ng-loading-bar',
  'ng-file-upload',
  'ui.bootstrap',
  'ui.router',
  'tooltip',
  'rzModule',
  './modules/charts/index',
  './modules/user/feedback-ctrl',
  './modules/user/help-ctrl',
  './modules/user/index',
  './modules/analysis/index',
  './modules/admin/index',
  './modules/common/util-service',
  './modules/common/project-service',
  './modules/common/poller-service',
  './modules/common/icon-directive',
  './modules/common/form-input-validate-directive',
  './modules/common/local-storage-polyfill',
  './modules/calibration/index',
  './modules/project/index',
  './modules/geospatial/index',
  './modules/programs/index',
  './modules/user/user-manager-service',
  './modules/ui/modal-service',
  './modules/ui/index'
], function (angular) {

  'use strict';

  return angular
    .module('app',
      [
        'angularFileUpload',
        'angular-loading-bar',
        'ui.bootstrap',
        'ui.router',
        'tooltip.module',
        'rzModule',
        'app.common.util-service',
        'app.common.form-input-validate',
        'app.common.project-service',
        'app.common.poller-service',
        'app.common.icon-directive',
        'app.feedback',
        'app.help',
        'app.user',
        'app.analysis',
        'app.admin',
        'app.charts',
        'app.local-storage',
        'app.model',
        'app.programs',
        'app.project',
        'app.geospatial',
        'app.ui',
        'app.ui.modal',
        'app.user-manager'
      ])

    .config(function ($httpProvider) {
      $httpProvider.interceptors.push(function ($q, $injector) {
        return {
          responseError: function (rejection) {
            if (rejection.status === 401 && rejection.config.url !== '/api/users/current') {
              // Redirect them back to login page
              location.href = './#/login';

              return $q.reject(rejection);
            } else {
              var message, errorText;
              console.log('catching error', rejection);
              if (rejection.data && (rejection.data.message || rejection.data.exception || rejection.data.reason)) {
                errorText = rejection.data.message || rejection.data.exception || rejection.data.reason;
              } else {
                errorText = JSON.stringify(rejection, null, 2);
              }
              message = 'We are very sorry, but it seems an error has occurred. Please contact us (info@optimamodel.com). In your email, copy and paste the error message below, and please also provide the date and time, your user name, the project you were working on (if applicable), and as much detail as possible about the steps leading up to the error. We apologize for the inconvenience.';
              var modalService = $injector.get('modalService');
              modalService.inform(angular.noop, 'Okay', message, 'Server Error', errorText);

              return $q.reject(rejection);
            }
          }
        };
      });
    })

    .config(function ($urlRouterProvider) {
      $urlRouterProvider.otherwise('/');
    })

    .run(function ($rootScope, $state, userManager, projectService) {

      /**
       * an injector has been run in main.js before app.js to fetch
       * the current user and stored in window.user, this will be
       * used to build the app in the first run
       */
      if (window.user) {
        userManager.setUser(window.user);
        delete window.user;
      }

      // Set the active project if any
      projectService.loadActiveProject();

      function isStatePublic(stateName) {
        var publicStates = ['contact', 'login', 'register'];
        return publicStates.indexOf(stateName) !== -1;
      }

      $rootScope.$on('$stateChangeStart', function (event, to) {
        if (!userManager.isLoggedIn && !isStatePublic(to.name)) {
          event.preventDefault();
          $state.go('login');
        }
      });
    });

});
