import React from 'react';
import {useNavigate} from "react-router-dom";

export default function ToolDescription({URL, imagePath, toolTitle, textDescription}) {

    const navigate = useNavigate();
    const navigateToToolPage = () => {
        navigate(URL);
    }

    return (
      <div className={'container tool-description'} onClick={navigateToToolPage}>
          <img className={'tool-reference-image'} src={imagePath} alt={'description'}/>
          <h3 className={'tool-title'}>{toolTitle}</h3>
          <p className={'tool-text-description'}>{textDescription}</p>
      </div>
    );
}