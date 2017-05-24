angular.module("demo").controller("NestedListsDemoController", function($scope) {

    $scope.models = {
        selected: null,
        templates: [
            {type: "item", typeName: "Component", id: 1},
            {type: "container", typeName: "Design", id: 2, columns: [[]]}
        ],
        designs: [
            {
                "type": "container",
                "id": 1,
                "columns": [[]]
            },
        ]
    };

    $scope.$watch('models.designs', function(model) {
        $scope.modelAsJson = angular.toJson(model, true);
    }, true);
});