designApp.directive("designPanel", function() {
    return {
    	scope: {
    		"models": "="
    	},
        templateUrl: "design.html"
    };
});