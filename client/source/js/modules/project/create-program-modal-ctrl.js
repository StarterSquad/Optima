define(['./module', 'angular', 'underscore'], function (module, angular, _) {
  'use strict';

  module.controller('ProjectCreateProgramModalController', function ($scope, $modalInstance, program) {

    // Initializes relevant attributes
    var initialize = function() {
      $scope.isNew = !program.name;
      $scope.program = program;
    };

    $scope.submit = function (form) {
      if (form.$invalid) {
        alert('Please fill in the form correctly');
      }

      $modalInstance.close(program);
    };

    initialize();

  });

});
