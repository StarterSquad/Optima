define(['angular', '../../version', 'ui.router'], function (angular, version) {
  'use strict';

  return angular.module('app.help', ['ui.router'])
    .config(function ($stateProvider) {
      $stateProvider
        .state('help', {
          url: '/help',
          templateUrl: 'js/modules/user/help.html',
          controller: function ($scope) {
            $scope.version = version;
          }
        });
    });

});
