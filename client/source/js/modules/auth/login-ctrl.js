define(['./module'], function (module) {
  'use strict';

  return module.controller('LoginController', function ($scope, $window, User) {

    $scope.error = '';

    $scope.login = function () {
      $scope.$broadcast('form-input-check-validity');

      if ($scope.LogInForm.$invalid) {
        return;
      }

      $scope.error = '';

      User.login({
        email: $scope.email,
        password: $scope.password
      },
        // success
        function (user) {
          $window.location = '/';
        },
        // error
        function (error) {
          if (error.status === 401) {
            $scope.error = 'Wrong email or password. Please check credentials and try again';
          } else {
            $scope.error = 'Server feels bad. Please try again in a bit';
          }
        }
      );
    };

  });

});
