designApp.directive("designPanel", function() {
    return {
    	scope: {
    		"addDesign": "&",
    		"models": "="
    	},
        templateUrl: "design.html"
    };
});