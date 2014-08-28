define([
    'angular',
    './button-choicebox/index',
    './menu/index'
], function (angular) {
    'use strict';

    return angular.module('app.ui', [
        'app.ui.button-choicebox',
        'app.ui.menu'
    ]).controller('MainCtrl', function ($scope) {
        $scope.asideMenuSettings = {
            items: [
                {
                    title: 'Home',
                    state: {
                        name: 'home'
                    }
                },
                {
                    title: 'Create/load project',
                    id: 'create-load',
                    subitems: [
                        {
                            title: 'Create new project',
                            state: {
                                name: 'project.create'
                            }
                        },
                        {
                            title: 'Upload data spreadsheet',
                            state: {
                                name: 'project.upload-data'
                            }
                        },
                        {
                            title: 'Load existing project',
                            state: {
                                name: 'project.load'
                            }
                        },
                        {
                            title: 'Modify existing project',
                            state: {
                                name: 'project.modify'
                            }
                        }
                    ]
                },
                {
                    title: 'View & calibrate model',
                    id: 'create-load',
                    subitems: [
                        {
                            title: 'View model & data',
                            state: {
                                name: 'model.view'
                            }
                        },
                        {
                            title: 'Automatic calibration',
                            state: {
                                name: 'model.automatic-calibration'
                            }
                        },
                        {
                            title: 'Manual calibration',
                            state: {
                                name: 'model.manual-calibration'
                            }
                        }
                    ]
                },
                {
                    title: 'Optimization',
                    id: 'create-load',
                    subitems: [
                        {
                            title: 'Define objectives',
                            state: {
                                name: 'optimization.objectives'
                            }
                        },
                        {
                            title: 'Define constraints',
                            state: {
                                name: 'optimization.constraints'
                            }
                        }
                    ]
                },
                {
                    title: 'View/output results',
                    state: {
                        name: 'results'
                    }
                }
            ]
        };

    });
});
