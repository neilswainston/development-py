designApp.controller("designCtrl", function($scope) {

    $scope.models = {
        selected: null,
        templates: [
            {type: "item", typeName: "Feature", id: 1},
            {type: "container", typeName: "Design", id: 2, "features": []}
        ],
        designs: [
            {
                "type": "container",
                "id": 1,
                "features": []
            },
        ]
    };

    $scope.$watch('models.designs', function(model) {
        $scope.modelAsJson = angular.toJson(model, true);
    }, true);
});