// Improved section highlighting for better accuracy
document.addEventListener('DOMContentLoaded', function() {
  // Get all section containers with IDs
  const sections = document.querySelectorAll('.section-container');
  
  // Get all navigation links that point to section IDs
  const navLinks = document.querySelectorAll('.navbar-nav .nav-link[href*="#"]');
  const homeLink = document.querySelector('.navbar-nav .nav-link[href="/"]');
  
  // Function to determine which section is most visible in the viewport
  function getMostVisibleSection() {
    let maxVisibleSection = null;
    let maxVisibleAmount = 0;
    
    sections.forEach(section => {
      const rect = section.getBoundingClientRect();
      const viewportHeight = window.innerHeight;
      
      // Calculate how much of the section is visible in the viewport
      const visibleTop = Math.max(0, rect.top);
      const visibleBottom = Math.min(viewportHeight, rect.bottom);
      const visibleHeight = Math.max(0, visibleBottom - visibleTop);
      
      // Factor in the position in the viewport (prefer sections at the top)
      const positionFactor = 1 - (visibleTop / viewportHeight) * 0.5;
      const visibleAmount = visibleHeight * positionFactor;
      
      if (visibleAmount > maxVisibleAmount) {
        maxVisibleAmount = visibleAmount;
        maxVisibleSection = section;
      }
    });
    
    return maxVisibleSection;
  }
  
  // Function to highlight the active section in the navigation
  function highlightActiveSection() {
    // Special case for top of page
    if (window.scrollY < 100) {
      navLinks.forEach(link => link.classList.remove('active'));
      if (homeLink) homeLink.classList.add('active');
      return;
    }
    
    // Get most visible section
    const activeSection = getMostVisibleSection();
    
    if (activeSection) {
      const sectionId = activeSection.getAttribute('id');
      
      // Remove active class from all links
      navLinks.forEach(link => link.classList.remove('active'));
      if (homeLink) homeLink.classList.remove('active');
      
      // Add active class to the corresponding link
      const activeLink = document.querySelector(`.nav-link[href*="#${sectionId}"]`);
      if (activeLink) {
        activeLink.classList.add('active');
      }
    }
  }
  
  // Handle Home link click
  if (homeLink) {
    homeLink.addEventListener('click', function(e) {
      e.preventDefault();
      
      window.scrollTo({
        top: 0,
        behavior: 'smooth'
      });
      
      history.pushState(null, null, '/');
      
      // Update active state for home
      navLinks.forEach(link => link.classList.remove('active'));
      homeLink.classList.add('active');
    });
  }
  
  // Add click event listeners to section links
  navLinks.forEach(link => {
    link.addEventListener('click', function(e) {
      const href = this.getAttribute('href');
      if (href.includes('#') && !href.startsWith('http')) {
        e.preventDefault();
        
        const targetId = href.split('#')[1];
        const targetSection = document.getElementById(targetId);
        
        if (targetSection) {
          // Set this link as active immediately
          navLinks.forEach(l => l.classList.remove('active'));
          if (homeLink) homeLink.classList.remove('active');
          this.classList.add('active');
          
          // Calculate position and scroll
          const navbarHeight = document.querySelector('#navbar').offsetHeight;
          const targetPosition = targetSection.getBoundingClientRect().top + window.pageYOffset - navbarHeight;
          
          window.scrollTo({
            top: targetPosition,
            behavior: 'smooth'
          });
          
          // After scrolling completes, verify correct active state
          setTimeout(() => {
            navLinks.forEach(l => l.classList.remove('active'));
            if (homeLink) homeLink.classList.remove('active');
            this.classList.add('active');
          }, 1000);
          
          history.pushState(null, null, href);
        }
      }
    });
  });
  
  // Listen for scroll events to update active section
  window.addEventListener('scroll', debounce(highlightActiveSection, 100));
  
  // Debounce function to limit the frequency of scroll event handling
  function debounce(func, wait) {
    let timeout;
    return function() {
      const context = this;
      const args = arguments;
      clearTimeout(timeout);
      timeout = setTimeout(() => func.apply(context, args), wait);
    };
  }
  
  // Initialize on page load
  highlightActiveSection();
});