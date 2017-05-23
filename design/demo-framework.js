angular.module("demo", ["ngRoute", "dndLists"])
    .config(function($routeProvider) {
        $routeProvider
            .when('/nested', {
                templateUrl: 'nested-frame.html',
                controller: 'NestedListsDemoController'
            })
            .otherwise({redirectTo: '/nested'});
    });