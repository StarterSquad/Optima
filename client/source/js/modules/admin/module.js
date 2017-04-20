define([
    'angular',
    'ui.router',
], function (angular) {
    'use strict';

    return angular.module('app.admin', [
        'ui.router'
    ]).config(function ($stateProvider) {
        $stateProvider
            .state('admin', {
                url: '/admin',
                abstract: true,
                template: '<div ui-view></div>'
            })
            .state('admin.manage-users', {
                url: '/manage-users',
                templateUrl: 'js/modules/admin/manage-users.html' ,
                controller: 'AdminManageUsersController',
                resolve: {
                  users:function($http){
                    return $http.get('/api/user');
                  }
                }
            })
            .state('admin.manage-projects', {
                url: '/manage-projects',
                templateUrl: 'js/modules/admin/manage-projects.html' ,
                controller: 'AdminManageProjectsController',
                resolve: {
                  projects: function (projectService) {
                    return projectService.getAllProjectList();
                  },
                  users:function($http){
                    return $http.get('/api/user');
                  }
                }
            });
    });
});
